function prd = computeShiftedScaledCopy(example, prm)


input = example * prm(1);

% shift the input
shift = round(prm(2)/0.001);
prd = [zeros(1, shift), input(1 : end - shift )']';

end