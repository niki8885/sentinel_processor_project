def pytest_configure(config):
    config.addinivalue_line("markers", "slow: mark test as slow")
    config.addinivalue_line("markers", "fortran: requires compiled Fortran library")
    config.addinivalue_line("markers", "integration: writes/reads real files")