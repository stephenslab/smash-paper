load('../input/ml_in.txt');
est = recbams(ml_in);
csvwrite('../input/ml_out.csv', est);
exit
