    trns = 1;   % i'th row to save in transformed response
    ts = 0.01; 
    time = [0:ts:100]';
    nt = length(time);
    fs = 1/ts;  % sampling frequency
    Nsamp = 10;  %change it to 10000 

    tic;
    dof = 5; 
    
%% Simulation

    mass = 10*eye(dof);
    stiff = [4000 -2000  0  0  0;
            -2000  4000 -2000  0  0;
             0 -2000  4000 -2000  0;
             0  0 -2000  4000 -2000;
             0  0  0 -2000  2000];

%% arrange according to the modes required, solving permutation ambiguity 
% ------------------------------------------------------------------------
    [evec, eval]= eig(mass\stiff);
    eta = 0.01; 
    Damp = inv(evec')*(2*eta*sqrt(eval))*inv(evec);
%     Damp= zeros(dof,dof);

    [evec] = permutation(dof,evec);

%% simulation starts here
% ------------------------------------------------------------------------
% evec= flipud(evec) 

for nsim = 1:1
    nsim    %% nsim = number of simulations                                                    

    accel = 2*randn(nt,1); %white noise

    N = length(accel);
    Y_dd = resp(dof,mass,stiff,Damp,accel,time); % Y_dd is acceleration response 

    resp_mat(:,:,nsim) = Y_dd; 
%% Response
    
end 
figure; 
subplot(511), plot(time,Y_dd(:,1));  title('First floor response'); 
subplot(512), plot(time,Y_dd(:,2));  title('Second floor response'); 
subplot(513), plot(time,Y_dd(:,3));  title('Third floor response'); 
subplot(514), plot(time,Y_dd(:,4));  title('Fourth floor response'); 
subplot(515), plot(time,Y_dd(:,5));  title('Roof response'); 

[f1 fabs1] = fftplot1(ts, Y_dd(:,1)); 
[f2 fabs2] = fftplot1(ts, Y_dd(:,2)); 
[f3 fabs3] = fftplot1(ts, Y_dd(:,3)); 
[f4 fabs4] = fftplot1(ts, Y_dd(:,4)); 
[f5 fabs5] = fftplot1(ts, Y_dd(:,5)); 

figure; 
subplot(511), plot(f1, fabs1);  title('First floor frequency spectrum');  xlim([0 4]); 
subplot(512), plot(f2, fabs2);  title('Second floor frequency spectrum');  xlim([0 4]); 
subplot(513), plot(f3, fabs3);  title('Third floor frequency spectrum');  xlim([0 4]); 
subplot(514), plot(f4, fabs4);  title('Fourth floor frequency spectrum');  xlim([0 4]); 
subplot(515), plot(f5, fabs5);  title('Roof frequency spectrum');  xlim([0 4]);

matx= Y_dd*evec; 
[f1x, fabs1x] = fftplot1(ts, matx(:,1)); 
[f2x, fabs2x] = fftplot1(ts, matx(:,2)); 
[f3x, fabs3x] = fftplot1(ts, matx(:,3)); 
[f4x, fabs4x] = fftplot1(ts, matx(:,4)); 
[f5x, fabs5x] = fftplot1(ts, matx(:,5)); 
figure; 
subplot(511), plot(f1x, fabs1x);  title('First floor frequency spectrum');  xlim([0 4]); 
subplot(512), plot(f2x, fabs2x);  title('Second floor frequency spectrum');  xlim([0 4]); 
subplot(513), plot(f3x, fabs3x);  title('Third floor frequency spectrum');  xlim([0 4]); 
subplot(514), plot(f4x, fabs4x);  title('Fourth floor frequency spectrum');  xlim([0 4]); 
subplot(515), plot(f5x, fabs5x);  title('Roof frequency spectrum'); xlim([0 4]);


%faulted data
%% Simulation

    mass = 10*eye(dof);
    stiff = [4000 -2000  0  0  0;
            -2000  4000 -2000  0  0;
             0 -2000  4000 -2000  0;
             0  0 -2000  4000 -2000;
             0  0  0 -2000  2000];
    k=0.7*stiff;

%% arrange according to the modes required, solving permutation ambiguity 
% ------------------------------------------------------------------------
    [evecf, evalf]= eig(mass\k);
    eta = 0.01; 
    Dampf = inv(evecf')*(2*eta*sqrt(evalf))*inv(evecf);
%     Damp= zeros(dof,dof);

    [evecf] = permutation(dof,evecf);

%% simulation starts here
% ------------------------------------------------------------------------
% evec= flipud(evec) 

for nsim = 1:1
    nsim    %% nsim = number of simulations                                                    

    accel = 2*randn(nt,1); %white noise

    N = length(accel);
    Y_ddf = resp(dof,mass,k,Dampf,accel,time); % Y_dd is acceleration response 

    resp_mat(:,:,nsim) = Y_ddf; 
%% Response
    
end 
figure; 
subplot(511), plot(time,Y_ddf(:,1));  title('First floor response'); 
subplot(512), plot(time,Y_ddf(:,2));  title('Second floor response'); 
subplot(513), plot(time,Y_ddf(:,3));  title('Third floor response'); 
subplot(514), plot(time,Y_ddf(:,4));  title('Fourth floor response'); 
subplot(515), plot(time,Y_ddf(:,5));  title('Roof response'); 

[f1 fabs1] = fftplot1(ts, Y_ddf(:,1)); 
[f2 fabs2] = fftplot1(ts, Y_ddf(:,2)); 
[f3 fabs3] = fftplot1(ts, Y_ddf(:,3)); 
[f4 fabs4] = fftplot1(ts, Y_ddf(:,4)); 
[f5 fabs5] = fftplot1(ts, Y_ddf(:,5)); 

figure; 
subplot(511), plot(f1, fabs1);  title('First floor frequency spectrum');  xlim([0 4]); 
subplot(512), plot(f2, fabs2);  title('Second floor frequency spectrum');  xlim([0 4]); 
subplot(513), plot(f3, fabs3);  title('Third floor frequency spectrum');  xlim([0 4]); 
subplot(514), plot(f4, fabs4);  title('Fourth floor frequency spectrum');  xlim([0 4]); 
subplot(515), plot(f5, fabs5);  title('Roof frequency spectrum');  xlim([0 4]);

matx= Y_ddf*evecf; 
[f1x fabs1x] = fftplot1(ts, matx(:,1)); 
[f2x fabs2x] = fftplot1(ts, matx(:,2)); 
[f3x fabs3x] = fftplot1(ts, matx(:,3)); 
[f4x fabs4x] = fftplot1(ts, matx(:,4)); 
[f5x fabs5x] = fftplot1(ts, matx(:,5)); 
figure; 
subplot(511), plot(f1x, fabs1x, 'r');  title('First floor frequency spectrum');  xlim([0 4]); 
subplot(512), plot(f2x, fabs2x, 'b');  title('Second floor frequency spectrum');  xlim([0 4]); 
subplot(513), plot(f3x, fabs3x, 'g');  title('Third floor frequency spectrum');  xlim([0 4]); 
subplot(514), plot(f4x, fabs4x, 'c');  title('Fourth floor frequency spectrum');  xlim([0 4]); 
subplot(515), plot(f5x, fabs5x, 'k');  title('Roof frequency spectrum'); xlim([0 4]);

%CorrCA start:


[m,n]= size(Y_dd);
cov_m= (Y_dd*matx')/n;
[U,S,V]=svd(cov_m);
figure;
plot(evec)
hold on
plot(U,'--')