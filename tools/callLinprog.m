function xsol = callLinprog(f,A,b,Aeq,beq,lb,ub,options)

v = version('-release');
yr = str2num(v(1:4));
if yr < 2016 || (yr == 2016 && v(end)=='a')
    xsol = linprog(f,A,b,Aeq,beq,lb,ub,[],options);
else
    xsol = linprog(f,A,b,Aeq,beq,lb,ub,options);
end
