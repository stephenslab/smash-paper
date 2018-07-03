load('../input/ml_in.txt');
est = recneighblock(ml_in);
csvwrite('../input/ml_out.csv', est);
exit


