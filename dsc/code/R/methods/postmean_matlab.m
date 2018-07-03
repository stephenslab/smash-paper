load('../input/ml_in.txt');
est = recsinglemean(ml_in, [0.1, 1]);
csvwrite('../input/ml_out.csv', est);
exit


