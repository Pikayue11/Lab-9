clear all;
global  initial_flag
initial_flag = 0;
D =10; d=D; %30,
Xmax = 100 + zeros(28,1);
Xmax(4) = 10;
Xmax(5) = 10;
Xmax(6) = 20;
Xmax(7) = 50;
Xmax(9) = 10;
Xmax(19) = 50;
Xmax(28) = 50;

Max_const_pro_FES=20000;

    tabulka=[];
    runs=25;
for func_num=1:28
    func_num
    a=-Xmax(func_num,1);
    b=Xmax(func_num,1);
    eps_viol=0.0001; 
     
    
    for j= 1:runs
        [vystup] = LSHADE44const(D,Max_const_pro_FES,a,b,func_num,eps_viol);
        vystup=[func_num j vystup];
        tabulka=[tabulka;vystup];     
    end
    
end
tabulka

    soubor=strcat('LSHADE44const_i_s_c_a_porusenim',num2str(D),'_behu',num2str(runs),'.txt');
    fid = fopen (soubor, 'wt');
    fprintf(fid,' %14.6g   %14.6g  %14.6g   %14.6g   %14.6g   %14.6g  %14.6g   %14.6g   %14.6g  %14.6g   %14.6g  %14.6g   %14.6g   %14.6g   %14.6g   %14.6g   %14.6g  %14.6g   %14.6g   %14.6g   %14.6g \n', tabulka');
    fclose(fid);

    
    clear all;
global  initial_flag
initial_flag = 0;
D =30; d=D; %30,
Xmax = 100 + zeros(28,1);
Xmax(4) = 10;
Xmax(5) = 10;
Xmax(6) = 20;
Xmax(7) = 50;
Xmax(9) = 10;
Xmax(19) = 50;
Xmax(28) = 50;

Max_const_pro_FES=20000;

    tabulka=[];
    runs=25;
for func_num=1:28
    func_num
    a=-Xmax(func_num,1);
    b=Xmax(func_num,1);
    eps_viol=0.0001; 
     
    
    for j= 1:runs
        [vystup] = LSHADE44const(D,Max_const_pro_FES,a,b,func_num,eps_viol);
        vystup=[func_num j vystup];
        tabulka=[tabulka;vystup];     
    end
    
end
tabulka

    soubor=strcat('LSHADE44const_i_s_c_a_porusenim',num2str(D),'_behu',num2str(runs),'.txt');
    fid = fopen (soubor, 'wt');
    fprintf(fid,' %14.6g   %14.6g  %14.6g   %14.6g   %14.6g   %14.6g  %14.6g   %14.6g   %14.6g  %14.6g   %14.6g  %14.6g   %14.6g   %14.6g   %14.6g   %14.6g   %14.6g  %14.6g   %14.6g   %14.6g   %14.6g \n', tabulka');
    fclose(fid);

    
    clear all;
global  initial_flag
initial_flag = 0;
D =50; d=D; %30,
Xmax = 100 + zeros(28,1);
Xmax(4) = 10;
Xmax(5) = 10;
Xmax(6) = 20;
Xmax(7) = 50;
Xmax(9) = 10;
Xmax(19) = 50;
Xmax(28) = 50;

Max_const_pro_FES=20000;

    tabulka=[];
    runs=25;
for func_num=1:28
    func_num
    a=-Xmax(func_num,1);
    b=Xmax(func_num,1);
    eps_viol=0.0001; 
     
    
    for j= 1:runs
        [vystup] = LSHADE44const(D,Max_const_pro_FES,a,b,func_num,eps_viol);
        vystup=[func_num j vystup];
        tabulka=[tabulka;vystup];     
    end
    
end
tabulka

    soubor=strcat('LSHADE44const_i_s_c_a_porusenim',num2str(D),'_behu',num2str(runs),'.txt');
    fid = fopen (soubor, 'wt');
    fprintf(fid,' %14.6g   %14.6g  %14.6g   %14.6g   %14.6g   %14.6g  %14.6g   %14.6g   %14.6g  %14.6g   %14.6g  %14.6g   %14.6g   %14.6g   %14.6g   %14.6g   %14.6g  %14.6g   %14.6g   %14.6g   %14.6g \n', tabulka');
    fclose(fid);

    
    clear all;
global  initial_flag
initial_flag = 0;
D =100; d=D; %30,
Xmax = 100 + zeros(28,1);
Xmax(4) = 10;
Xmax(5) = 10;
Xmax(6) = 20;
Xmax(7) = 50;
Xmax(9) = 10;
Xmax(19) = 50;
Xmax(28) = 50;

Max_const_pro_FES=20000;

    tabulka=[];
    runs=25;
for func_num=1:28
    func_num
    a=-Xmax(func_num,1);
    b=Xmax(func_num,1);
    eps_viol=0.0001; 
     
    
    for j= 1:runs
        [vystup] = LSHADE44const(D,Max_const_pro_FES,a,b,func_num,eps_viol);
        vystup=[func_num j vystup];
        tabulka=[tabulka;vystup];     
    end
    
end
tabulka

    soubor=strcat('LSHADE44const_i_s_c_a_porusenim',num2str(D),'_behu',num2str(runs),'.txt');
    fid = fopen (soubor, 'wt');
    fprintf(fid,' %14.6g   %14.6g  %14.6g   %14.6g   %14.6g   %14.6g  %14.6g   %14.6g   %14.6g  %14.6g   %14.6g  %14.6g   %14.6g   %14.6g   %14.6g   %14.6g   %14.6g  %14.6g   %14.6g   %14.6g   %14.6g \n', tabulka');
    fclose(fid);
