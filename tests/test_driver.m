addpath('..\src\');
suite = matlab.unittest.TestSuite.fromClass(?TestLoadDataset);
result = suite.run