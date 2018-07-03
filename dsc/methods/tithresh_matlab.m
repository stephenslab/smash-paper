load('../input/ml_in.txt');
est = recTI(ml_in);
csvwrite('../input/ml_out.csv',est);
exit
