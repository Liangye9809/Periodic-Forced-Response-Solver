% plot backwards from 0.9 on each dof
load('A_preload_backwards_cont_friction_F11000_1H_without_initialw.mat');
for i = 1:5
    figure(i)
    plot(A_cont(:,6), A_cont(:,i), 'b-'), hold on;
end

% plot forwards from 0.9 on each dof
load('A_preload_forward_cont_friction_F11000_1H_without_initialw.mat');
for i = 1:5
    figure(i)
    plot(A_cont(:,6), A_cont(:,i), 'r-'), hold on;
end


% plot backwards from 1.2 on each dof
load('A_cont_Backwards_H_10.mat');
for i = 1:5
    figure(i)
    plot(A_cont(:,6), A_cont(:,i), 'k--'), hold on;
end

% plot forwards from 1.2 on each dof
load('A_cont_Forwards_H_10.mat');
for i = 1:5
    figure(i)
    plot(A_cont(:,6), A_cont(:,i), 'k--'), hold on;
end

% legend of each plot
for i = 1:5
    figure(i)
    legend('backwards from 1.3','forwards from 0.9','backwards from 1.2','forwards from 1.2');
    grid on;
    xlabel('Omega');
    if i < 3
        name = 'q_' + string(i);
        ylabel('Amplitude of' + name);
        title(name + ' - Omega');
    else
        name = 'x_' + string(i-2);
        ylabel('Amplitude of' + name);
        title(name + ' - Omega');
    end
end

%%
% friction forces over time
cd data/'27 May - detail in omega 1.2'/
load("Ff_omega1.2_forwards.mat");
load("w_omega1.2_forwards.mat");
cd ../
cd ../
Ff_f = hndn;
Dissip_f = zeros(2,1);
% dissipation by Harmonics' coefficient
for i = 1:2
    Dissip_f(i,1) = 0.5 * T * 1.2 * (Xforwards(3:2:end, i + 2)' * Ff_f(2:2:end, i) -  Ff_f(3:2:end, i)' * Xforwards(2:2:end, i + 2));
end

ft_f = E * Ff_f - params.func.static.preload.gxp'; % g(x + xp)
w_f = w_cont;
figure;
for i = 1:2
    plot(t, ft_f(:,i)), hold on;
end
plot(t, ft_f(:,3), 'k-'), hold on;
plot(t, -ft_f(:,3), 'k-'), hold on;
legend('T1', 'T2', 'Fn');
grid on;
xlabel('t');
ylabel('Force');
title('Friction forces over time forwards');

cd data/'27 May - detail in omega 1.2'/
load("Ff_omega1.2_backwards.mat");
load("w_omega1.2_backwards.mat");
cd ../
cd ../
Ff_b = hndn;
Dissip_b = zeros(2,1);
% dissipation by Harmonics' coefficient
for i = 1:2
    Dissip_f(i,1) = 0.5 * T * 1.2 * (Xbackwards(3:2:end, i + 2)' * Ff_b(2:2:end, i) -  Ff_b(3:2:end, i)' * Xbackwards(2:2:end, i + 2));
end
ft_b = E * Ff_b - params.func.static.preload.gxp'; % g(x + xp)
w_b = w_cont;
figure;
for i = 1:2
    plot(t, ft_b(:,i)), hold on;
end
plot(t, ft_b(:,3),'k-'), hold on;
plot(t, -ft_b(:,3),'k-'), hold on;
legend('T1', 'T2', 'Fn');
grid on;
xlabel('t');
ylabel('Force');
title('Friction forces over time backwards');