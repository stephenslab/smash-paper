load('../input/ml_in.txt');
est = recsure(ml_in);
csvwrite('../input/ml_out.csv', est);
exit
