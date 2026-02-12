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

%% find transition time, from slip or saperation to stick
Mf = [1 0 1 1 1 0 -1 -1 0 1  1  1]';
ll = [1 2 3 4 5 6  7  8 9 10 11 12];
N = length(Mf);
record = zeros(size(Mf));
p = 0;
keep = 0;
for i = 1:N

    if i == 1 && Mf(1) ~=0 
        for j = 1:N
            k = N - j + 1;
            if Mf(k) == 0
                p = mod(k, N) + 1;
                keep = 1;
                break;
            end
        end
    elseif Mf(i) ~= 0 && keep == 0
        p = i;
        keep = 1;
    end

    if Mf(i) == 0 && keep == 1
        keep = 0;
    end

    record(i) = p;


end
record