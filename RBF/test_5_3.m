%%% 5.3节
% 径向基函数配置法求四种情形的问题

%% Omega_1 + u_1
clc
clear

bdyfunc = @bdyfunc1;
[idx, errmax, errmsq] = rbflaplace(bdyfunc,5);
figure(1)
semilogy(idx,errmax,'.-',idx,errmsq,'s-')
legend('errmax','errmsq')
xlabel('N')
ylabel('error')

%% Omega_1 + u_2
clc
clear

bdyfunc = @bdyfunc2;
[idx, errmax, errmsq] = rbflaplace(bdyfunc,5);
figure(2)
semilogy(idx,errmax,'.-',idx,errmsq,'s-')
legend('errmax','errmsq')
xlabel('N')
ylabel('error')

%% Omega_2 + u_1
clc
clear

bdyfunc = @bdyfunc1;
[idx, errmax, errmsq] = rbflaplace(bdyfunc,3);
figure(3)
semilogy(idx,errmax,'.-',idx,errmsq,'s-')
legend('errmax','errmsq')
xlabel('N')
ylabel('error')

%% Omega_2 + u_2
clc
clear

bdyfunc = @bdyfunc2;
[idx, errmax, errmsq] = rbflaplace(bdyfunc,3);
figure(4)
semilogy(idx,errmax,'.-',idx,errmsq,'s-')
legend('errmax','errmsq')
xlabel('N')
ylabel('error')

