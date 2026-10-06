""" Tests for reading ECSV files with several fragments: the GFE standard (Appendix 4) keeps one row per
    frame, with the additional fragments in columns with a numeric suffix (e.g. azimuth1). Frames on which
    only an additional fragment was measured leave the main fragment columns empty, and main fragment points
    flagged with trajectory_use = False are left out of its trajectory.

Run with pytest:
    python -m pytest wmpl/Formats/Tests/test_ECSV.py -v
"""

from __future__ import print_function, division, absolute_import

import numpy as np

from wmpl.Formats.ECSV import ecsvFragments, loadECSVs


META = """# %ECSV 0.9
# ---
# meta: !!omap
# - {obs_latitude: -43.6}
# - {obs_longitude: 172.1}
# - {obs_elevation: 100.0}
# - {camera_id: 'STATION'}
# - {image_file: 'FF_STATION.fits'}
"""

MAIN = ['datetime', 'ra', 'dec', 'azimuth', 'altitude', 'x_image', 'y_image', 'mag_data']
FRAGMENT = ['datetime', 'ra', 'dec', 'azimuth', 'altitude', 'x_image', 'y_image']


def point(i, altitude):
    """ Time and position of a point on frame i. """

    return ['2026-01-30T10:25:36.{:06d}'.format(40000*i), '10.0', '-30.0', '{:.4f}'.format(120 + 0.5*i), 
        '{:.4f}'.format(altitude - 0.3*i), '{:.3f}'.format(100 + 10*i), '{:.3f}'.format(200 + 5*i)]


def writeEcsv(path, station, fragments=True):
    """ Six frames of the main fragment, the third flagged as not to be used and the last only measured for
        fragment 1, and fragment 1 on the last five frames. Without fragments, the previous file format. """

    header = MAIN + (['frame_number', 'flare', 'trajectory_use'] + [name + '1' for name in FRAGMENT] 
        if fragments else [])

    rows = []
    for i in range(7):

        main = point(i, 30.0) + ['{:.2f}'.format(2.0 - 0.1*i)] if i < 6 else [''] * len(MAIN)
        if not fragments:
            if i < 6:
                rows.append(main)
            continue

        flags = [str(i), 'False', 'False' if i == 2 else 'True'] if i < 6 else [str(i), '', '']
        fragment = point(i, 29.0) if i >= 2 else [''] * len(FRAGMENT)
        rows.append(main + flags + fragment)

    with open(path, 'w') as f:
        f.write(META.replace('STATION', station) + ",".join(header) + "\n")
        f.write("".join(",".join(row) + "\n" for row in rows))

    return str(path)


def test_fragments_are_found(tmp_path):
    """ The additional fragments are the suffixes of the azimuth columns; files without them have none. """

    assert ecsvFragments([writeEcsv(tmp_path/"a.ecsv", "XX0001")]) == [1]
    assert ecsvFragments([writeEcsv(tmp_path/"b.ecsv", "XX0002", fragments=False)]) == []


def test_main_fragment_skips_empty_and_unused_points(tmp_path):
    """ Of the 7 rows, the frame with only fragment 1 is skipped, and so is the main fragment point flagged
        as not to be used in the trajectory, leaving 5 points. """

    _, meteors = loadECSVs([writeEcsv(tmp_path/"a.ecsv", "XX0001")], no_prepare=True)

    assert len(meteors) == 1
    assert np.allclose(np.degrees(meteors[0].azim_data), [120.0, 120.5, 121.5, 122.0, 122.5])


def test_additional_fragment_columns(tmp_path):
    """ Fragment 1 is read from its own columns, without the photometry of the main fragment. """

    _, meteors = loadECSVs([writeEcsv(tmp_path/"a.ecsv", "XX0001")], no_prepare=True, fragment=1)

    assert len(meteors) == 1
    assert len(meteors[0].time_data) == 5
    assert np.isclose(np.degrees(min(meteors[0].elev_data)), 29.0 - 0.3*6)
    assert set(meteors[0].mag_data) == {10.0}


def test_files_without_the_fragment_are_skipped(tmp_path):
    """ Only the stations whose files describe a fragment are loaded for it. """

    paths = [writeEcsv(tmp_path/"a.ecsv", "XX0001"), writeEcsv(tmp_path/"b.ecsv", "XX0002", fragments=False)]

    _, main = loadECSVs(paths, no_prepare=True)
    _, fragment = loadECSVs(paths, no_prepare=True, fragment=1)

    assert sorted(meteor.station_id for meteor in main) == ['XX0001', 'XX0002']
    assert [meteor.station_id for meteor in fragment] == ['XX0001']


def test_previous_format_is_unchanged(tmp_path):
    """ A file without the annotation columns loads all its points, as before. """

    _, meteors = loadECSVs([writeEcsv(tmp_path/"b.ecsv", "XX0002", fragments=False)], no_prepare=True)

    assert len(meteors[0].time_data) == 6
    assert list(meteors[0].mag_data) == [2.0, 1.9, 1.8, 1.7, 1.6, 1.5]
