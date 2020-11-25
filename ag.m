close all;
clear;
clc;
%-------------------------------------
%Podanie wartości początkowych z przedziału (-10,10)
przedzial = (-10:0.1:10);
disp('Podaj wartośi z przedziału (-10,10)')
a0 = -11;
while a0 <- 10 || a0 > 10
    a0 = input('Podaj wartośc a0:');
end
a1 = -11;
while a1 <- 10 || a1 > 10
    a1 = input('Podaj wartośc a1:');
end
a2 = -11;
while a2 <- 10 || a2 > 10
    a2 = input('Podaj wartośc a2:');
end
b0 = -11;
while b0 <- 10 || b0 > 10
    b0 = input('Podaj wartośc b0:');
end
b1 = -11;
while b1 <- 10 || b1 > 10
    b1 = input('Podaj wartośc b1:');
end
%Rozmiar populacji
population_size = 2000;
%Ilość osobników pobrana do crrosowania
population_pick_cross = population_size / 2;
%ilosc generacji
generations = 20;
%Szansa na inwersje przy crossowaniu
inversion_rate = 0.01;
%mutacja
mutation_rate = 0.1;
%Ilość uruchomień programu
execute=5;
%---------------
%Przedział (0,100)
om = logspace(-2, 2, 401);
%Stworzenie transmitancji na podstawie danych początkowych podanych przez
%użytkownika
num = [a2 a1 a0];
den = [1 b1 b0];
sys = tf(num, den)
G = freqresp(sys, om);
G = G(:);
bodeG = bode(sys, om);
nyquistG = nyquist(sys, om);
%%
genes = generate_rand(); %gotowe
tablica_struct = [];
tablica_fitness = [];
tablica_error = [];
tablica_best_score = [];
tablica_error_best = [];
tablica_new_input=zeros(5,5);
%Stworzenie początkowej populacji
for l=1:execute
   %Pomocnicza zmienna do wyświetlania aktualnej generacji
   generation_num=0; 
for i = 1:population_size
    s.genes = generate_rand();
    tablica_struct = [tablica_struct s];
end


%Zapis fitnessów do tablicy
tablica_fitness = get_fitness_vector(tablica_struct, om, G, nyquistG, bodeG);

for p = 1:generations

    %Wypisanie aktualnej generacji
    
    generation_num = generation_num + 1

    %Ustawienie wartości od najmniejszej do największej
    [fitness_sorted, order] = sort(tablica_fitness);
    struct_sorted = tablica_struct(order);

    %Wybranie połowy populacji z lepszym fitnessem
    half_fitness_sorted = fitness_sorted(population_pick_cross:population_size);
    half_struct_sorted = struct_sorted(population_pick_cross:population_size);

    new_guys = [];

    %Crossowanie populacji
    for o = 1:population_pick_cross
        %wybór rodziców do crossowania
        parents = randsample(half_struct_sorted, 2, true, half_fitness_sorted);
        %Stworzenie nowego osobnika
        new_guy = get_new_guy(parents(1), parents(2), inversion_rate);
        %Zapis nowych osobników do tablicy
        new_guys = [new_guys new_guy];
    end
    %Dodanie mutacji do populacji, zapobiega to zatrzymaniu się algorytmu 
    for o=1 : length(new_guys) 
        new_guys(o) = mutation_flip(new_guys(o),mutation_rate);
    end
    
    %Obliczenie fitnessu dla nowych osobników
    new_fitness = get_fitness_vector(new_guys, om, G, nyquistG, bodeG);

    %Dopisanie do połowy starych osobników drugiej połowy nowych
    tablica_struct = [half_struct_sorted new_guys];
    %Dopisanie do połowy fitnessów starych osobników drugiej połowy nowych
    %fitnessów
    tablica_fitness = [half_fitness_sorted new_fitness];
    
    
    half_fitness_sorted(1)
end

%najlepszy osobnik
best_guy = tablica_struct(1);

%Obliczenie nowych wartości do transmitancji (a0,a1,a2,b0,b1)
%Dla jednego współczynnika przypada tablica [1,11]
new_a0 = best_guy.genes(1, :);
new_a0 = translate_coeff(-10, 10, 11, new_a0);

new_a1 = best_guy.genes(2, :);
new_a1 = translate_coeff(-10, 10, 11, new_a1);

new_a2 = best_guy.genes(3, :);
new_a2 = translate_coeff(-10, 10, 11, new_a2);

new_b0 = best_guy.genes(4, :);
new_b0 = translate_coeff(-10, 10, 11, new_b0);

new_b1 = best_guy.genes(5, :);
new_b1 = translate_coeff(-10, 10, 11, new_b1);

%Stworzenie transmitancji dla najlepszego osobnika
new_best_num = [new_a2 new_a1 new_a0];
new_best_den = [1 new_b1 new_b0];
new_best_sys = tf(new_best_num, new_best_den);

new_best_G = freqresp(new_best_sys, om);
new_best_G = new_best_G(:);
tablica_best_score(l) = tablica_fitness(1);
tablica_error_best(l)=1/tablica_fitness(1);
tablica_new_input(:,l) = [new_a0 ; new_a1 ; new_a2; new_b0;  new_b1];

% new_best_G_bode = bode(new_best,om);
% new_best_G_nyquist = nyquist(new_best,om);

%Wykresy porównujące nową transmitancję z początkową
if l==1
    %nyquist
figure(1)
plot(G,'b-','LineWidth',2)
hold on;
plot(new_best_G,'--','LineWidth',1)
hold on;

%bode
figure(2)
bode(sys)
hold on;
bodemag(new_best_sys)
hold on;

%skokowa
figure(3)
step(sys)
hold on;
step(new_best_sys)
hold on;

%best fitness
figure(4)
hold on;
plot(tablica_best_score,'r.','MarkerSize',20)
hold on;

%best a0,a1,a2,b0,b1
figure(5)
plot(tablica_new_input,'b.','MarkerSize',20)


else 
figure(1)
hold on;
plot(new_best_G,'--','LineWidth',1)
hold on;

figure(2)
hold on;
bodemag(new_best_sys)
hold on;

figure(3)
hold on;
step(new_best_sys)
hold on;

figure(4)
hold on;
plot(tablica_best_score,'r.','MarkerSize',20)
hold on;

figure(5)
plot(tablica_new_input,'b.','MarkerSize',20)

end
if l==execute
    figure(1)
    title('Charakterystyka częstotliwościowa')
    xlabel('Re')
    ylabel('Im')
    legend('Transmitancja zadana','data1','data2','data3','data4','data5')
    
    figure(2)
    legend('Transmitancja zadana','data1','data2','data3','data4','data5')
    
    figure(3)
    legend('Transmitancja zadana','data1','data2','data3','data4','data5')
    
    figure(4)
    title('Best Score')
    xlabel('Numer uruchomienia programu')
    ylabel('Najlepsza uzyskana wartość fitness')
    legend('Best Score')
    
    figure(5)
    axis([0 execute+1 -10 10])
    title('Wartości a0,a1,a2,b0,b1 dla najlepszego osobnika')
end
    
end
%Obliczanie mediany, wartości średniej oraz wariancji uzyskanych wyników 
median_error = median(tablica_error_best(1,:));

disp('Mediana dla błędu = ')
disp(median_error);


avarage_error = mean(tablica_error_best(1,:));

disp('Wartość średnia błędu = ')
disp(avarage_error);

variance_error = var(tablica_error_best(1,:));


disp('Wariancja błędu = ')
disp(variance_error);


%Funkcja wprowadzająca mutację do populacji
function mutation_flip = mutation_flip(guy,mutation_chance)


r= rand;
   if r < mutation_chance      
    genes = reshape(guy.genes,1,[]);
    genes_length=length(genes);
    a = randi(genes_length - 1);
    b = randi(genes_length - a)+a;
    fragment = genes(1,a:b);
    genes(1,a:b) = flip(fragment);
    guy.genes = reshape(genes,5,[]);
   end
mutation_flip = guy;

end
%Funkcja odpowiedzialna za obliczenie fitnessów
function get_fitness_vector = get_fitness_vector(guys, om, G, nyquistG, bodeG)

    tablica_fitness = [];

    for j = 1:length(guys)
        s = guys(j);
        a0 = translate_coeff(-10, 10, 11, s.genes(1, :));
        a1 = translate_coeff(-10, 10, 11, s.genes(2, :));
        a2 = translate_coeff(-10, 10, 11, s.genes(3, :));
        b0 = translate_coeff(-10, 10, 11, s.genes(4, :));
        b1 = translate_coeff(-10, 10, 11, s.genes(5, :));

        new_transmitancja = transmitancja(a0, a1, a2, b0, b1);
        G_new = freqresp(new_transmitancja, om);
        G_new = G_new(:);
        bodeG_new = bode(new_transmitancja, om);
        nyquistG_new = nyquist(new_transmitancja, om);

        g_error_freq = norm(G_new - G, 2);
        g_error_nyq = norm(nyquistG_new(1, :) - nyquistG(1, :), 2);
        g_error_bode = norm(bodeG_new(1, :) - bodeG(1, :), 2);
        
%Obliczanie błędu na podstawie uzyskanych wartości 
        error_new = (g_error_freq) + (g_error_bode) + (g_error_nyq);
        
%Obliczanie fitnessu
        fitness = 1 / error_new;

        tablica_fitness(j) = fitness;
    end

    get_fitness_vector = tablica_fitness;
end

%Funckja odpowiedzialna za stworzenie nowego osobnika
function get_new_guy = get_new_guy(parent1, parent2, inversion_rate)

    %wpisanie do nowego osobnika genotypu rodzica1
    new_guy = parent1;
    %zastąpienie 3 wierszy genotypu nowego osobnika 3 wierszami rodzica2
    new_guy.genes(3:5, :) = parent2.genes(3:5, :);

    for k = 1:5

        for l = 1:11
            a = rand;

            if a < inversion_rate
                new_guy.genes(k, l) = ~new_guy.genes(k, l);
            end

        end

    end

    get_new_guy = new_guy;
end

%Funkcja odpowiedzialna za stworzenie nowej transmitancji
function transmitancja = transmitancja(a0, a1, a2, b0, b1)

    num_new = [a2 a1 a0];
    den_new = [1 b1 b0];
    transmitancja = tf(num_new, den_new);

end

%Funckja tłumacząca fragment genotypu na współczynnik
function translate_coeff = translate_coeff(a, b, L, c)
    ssum = 0;

    for i = 0:L - 1
        ssum = ssum + c(:, i + 1) * 2^i;
    end

    x = a + (b - a) / (2^L - 1) * ssum;
    translate_coeff = x;
end


%% DONE
%Funckja generująca losowy genotyp
function generate_rand = generate_rand()

    u1 = randi([0, 1], 5, 11);

    generate_rand = u1;

end
