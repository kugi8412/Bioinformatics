% Parametry początkowe
K1 = 1;
K2 = 2;
N1_0 = 0.8;
N2_0 = 1.2;
a12_values = linspace(0.25, 1, 25);
a21_values = linspace(0.5, 1.5, 25);
threshold = 0.001;

% Zdarzenie wymarcia
function [value, isterminal, direction] = extinction_event(t, N)
    value = min(N) - (10 * eps);
    isterminal = 1;
    direction = 0;
end

% Siatka parametrów
grid_a12 = [];
grid_a21 = [];
survival_colors = [];

for a12 = a12_values
    for a21 = a21_values

        % Funkcja różniczkowa
        diff = @(t, N) [N(1)*(1 - N(1)/K1) - a12*N(1)*N(2);
                              N(2)*(1 - N(2)/K2) - a21*N(1)*N(2)];
        options = odeset('RelTol', 1e-9, 'AbsTol', 1e-9, 'MaxStep', 1, ...
                         'Events', @extinction_event);

        [t, N, te, ye, ie] = ode45(diff, [0 1000], [N1_0; N2_0], options);

        % Sprawdzenie wartości populacji
        if isempty(te)
            N1_final = N(end, 1);
            N2_final = N(end, 2);
        else
            N1_final = ye(end, 1);
            N2_final = ye(end, 2);
        endif

        % Przetrwanie populacji
        if N1_final > threshold  && N2_final < threshold
            survival_colors = [survival_colors; 1];
        elseif N1_final < threshold && N2_final > threshold
            survival_colors = [survival_colors; 2];
        elseif N1_final > threshold && N2_final > threshold
            survival_colors = [survival_colors; 3];
        else
            survival_colors = [survival_colors; 0];
        endif

        grid_a12 = [grid_a12; a12];
        grid_a21 = [grid_a21; a21];
    endfor
endfor

% Rysowanie wykresu
figure(1);
hold on;

% Wymuszenie pełnej legendy
plot(NaN, NaN, 'ro', 'MarkerSize', 5);
plot(NaN, NaN, 'bo', 'MarkerSize', 5);
plot(NaN, NaN, 'go', 'MarkerSize', 5);
plot(NaN, NaN, 'ko', 'MarkerSize', 5);

for i = 1:length(grid_a12)
    if survival_colors(i) == 1
        plot(grid_a12(i), grid_a21(i), 'ro', 'MarkerSize', 2);
    elseif survival_colors(i) == 2
        plot(grid_a12(i), grid_a21(i), 'bo', 'MarkerSize', 2);
    elseif survival_colors(i) == 3
        plot(grid_a12(i), grid_a21(i), 'go', 'MarkerSize', 2);
    else
        plot(grid_a12(i), grid_a21(i), 'ko', 'MarkerSize', 2);
    endif
endfor

xlabel('a_{12}');
ylabel('a_{21}');
title('Przetrwanie populacji');
legend('N_1 przetrwała', 'N_2 przetrwała', 'Obie przetrwały', ...
       'Żadna nie przetrwała','Location', 'bestoutside')
hold off;

% Wykresy fazowe dla przykładowych parametrów
example_params = [0.4, 0.8; 0.9, 1.2; 0.3, 1.4];

for i = 1:size(example_params, 1)
    a12 = example_params(i, 1);
    a21 = example_params(i, 2);

    % Funkcja różniczkowa
    ode_system = @(t, N) [N(1)*(1 - N(1)/K1) - a12*N(1)*N(2);
                          N(2)*(1 - N(2)/K2) - a21*N(1)*N(2)];

    % Wybrane rozwiązania
    options = odeset('RelTol', 1e-9, 'AbsTol', 1e-9, 'MaxStep', 0.1, ...
                     'Events', @extinction_event);

    [t, N, te, ye, ie] = ode45(ode_system, [0 2000], [N1_0; N2_0], options);

    % Wykres fazowy
    figure(i+1);
    plot(N(:,1), N(:,2), 'LineWidth', 1.5);
    hold on;
    plot(N(1,1), N(1,2), 'ro', 'MarkerSize', 5);
    if ~isempty(te)
        plot(ye(1), ye(2), 'ks', 'MarkerSize', 5);
    else
        plot(N(end,1), N(end,2), 'ks', 'MarkerSize', 5);
    endif

    xlabel('N_1');
    ylabel('N_2');
    title(sprintf('Wykres fazowy dla a_{12} = %.2f, a_{21} = %.2f', a12, a21));
    legend('Populacja', 'Początek', 'Koniec');
    hold off;
endfor

