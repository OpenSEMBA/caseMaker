import os

TEST_DATA_FOLDER = os.path.join(os.getcwd(), 'testData/')
CASES_FOLDER = os.path.join(TEST_DATA_FOLDER, 'cases/')
EXCITATIONS_FOLDER = os.path.join(TEST_DATA_FOLDER, 'excitations/')
OUTPUTS_FOLDER = os.path.join(TEST_DATA_FOLDER, 'outputs/')
GEOMETRIES_FOLDER = os.path.join(TEST_DATA_FOLDER, 'geometries/')
THIN_SLOT_FOLDER = os.path.join(TEST_DATA_FOLDER, 'thin_slot/')
THIN_SLOT_SURFACE_FOLDER = os.path.join(THIN_SLOT_FOLDER, 'surface/')


def find_fdtd_executable():
    repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    candidates = [
        os.environ.get('SEMBA_EXE'),
        os.path.join(repo_root, 'fdtd/build/bin/semba-fdtd'),
        os.path.join(repo_root, 'fdtd/build-rls/bin/semba-fdtd'),
        os.path.join(repo_root, 'fdtd/build-dbg/bin/semba-fdtd'),
    ]
    for candidate in candidates:
        if candidate and os.path.isfile(candidate) and os.access(candidate, os.X_OK):
            return candidate
    return None
