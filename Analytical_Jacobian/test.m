skip = [1 2; 2 1; 3 1; 3 2];
for i = 1:3
    for j = 1:3
        if ismember([i j], skip, 'rows')
            continue;
        end
        disp([i,j]);
    end
end
%% plot 2D tangential relative displacement

T = 1;
N = 256;
phy = pi / 2;
dt = T / N;
t = 0:dt:10*pi;

figure;
plot(sin(t),cos(t), 'k-'), hold on;
grid on;
for Uy0 = 0.1:0.2:1
% Uy0 = 1;
    Ux0 = 2.5 * Uy0;
    
    Ux = Ux0 * sin(2*pi*t);
    Uy = Uy0 * sin(2*pi*t + phy);
    % figure;
    plot(Ux, Uy), hold on;

% drawnow;
% pause(0.05);
end
xlim([-3,3]);
ylim([-3,3]);
pbaspect([1 1 1]);

%% find transition time, from slip or separation to stick
% a = [1 1 0 0 0 0 0 0 1 1 1 1 1 1 0 -1 -1 0 1  1  1]'; % 0 is slip or separation. 1 is stick
a(:, 1) = Mft_A(1,1,:);
b(:, 1) = Mft_A(2,1,:);
c(:, 1) = Mft_A(3,1,:);
ll = 1:size(a,1);
N = length(a);
record = zeros(size(a));
p = 0;
keep = 0;
for i = 1:N

    if i == 1 && a(1) ~=0 
        for j = 1:N
            k = N - j + 1;
            if a(k) == 0
                p = mod(k, N) + 1;
                keep = 1;
                break;
            end
        end
    elseif a(i) ~= 0 && keep == 0
        p = i;
        keep = 1;
    end

    if a(i) == 0 && keep == 1
        keep = 0;
        p = 0;
    end

    record(i) = p;


end
c_xn = a .* b(mod(record-2, N));

conditi = [a, c_xn];