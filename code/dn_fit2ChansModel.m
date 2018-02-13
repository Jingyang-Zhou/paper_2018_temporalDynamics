function target = dn_fit2ChansModel(params, data, t, stimulus, dt, fields)


x = [];
x      = toSetField(x, fields, params);
output = dn_2Chansmodel(x, stimulus, t, dt);


target = sum((output - data).^2);

%figure (90), clf, plot(data, 'r-'), hold on, plot(output, 'b-'), drawnow
end