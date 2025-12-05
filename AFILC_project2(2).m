clc,clear,close all;

% 系統參數
m = 2.3;          % 質量 (kg)
l = 0.9;          % 長度 (m)
J = 1.1667;       % 慣量矩 (kg·m^2)
g = 9.8;          % 重力加速度 (m/s^2)
T = 12;           % 模擬總時間 (s)
Td = 0.4;         % 時間步長
K = 10;           % 最大迭代次數
c_1 = 1;
lambda_1 = 6;
gamma_1 = 0.9;
gamma_2 = 0.6;
beta_1 = 2.5;
beta_2 = 12;
kc_upper = 0.2;
kc_lower = 0.15;
sp = 15;
dp = 1.5;
sn = -15;
dn = -1.5;
u_max = 15;
u_min = -15;
L = [6,5,4,3,2,1,0,-1,-2,-3,-4,-5,-6]';

% 控制器參數
lambda1 = 6;     % 控制增益
gamma1 = 0.9;    
beta1 = 2.5;
gamma2 = 0.6;
beta2 = 2.0;

% 時間設定
dt = 0.003;
time = 0:dt:T;

% 儲存狀態變量
x = zeros(2, length(time)+1); 

% 定義期望軌跡 (目標)
xd_1 = @(j,t) ((0.3 + 0.01*(j/K)).*sin((1 + 0.01*(j/K))* t)) + 0.2*cos(t); % 期望的角度
xd_2 = @(j,t) ((0.3+0.01*(j/K)).*(1+0.01*(j/K)).*cos((1+0.01*(j/K))*t)) - 0.2*sin(t); % 期望的速度
xd_3 = @(j,t) -(((0.3+0.01*(j/K)).*(1+0.01*(j/K))^2).*(sin((1+0.01*(j/K))*t)))-0.2*cos(t);

% 模糊隸屬度函數
L_centers = [6;5;4;3;2;1;0;-1;-2;-3;-4;-5;-6];
psi_fls_M1 = @(X,L) exp(-(X - L).^2/10)';

% 儲存控制輸入
u = zeros(1, length(time));
x1_store = zeros(1, length(time), 10);
x2_store = zeros(1, length(time), 10);
w_store = zeros(length(time),13,10);
xd_store = zeros(3,length(time),10);
w_hat_diff_store = zeros(length(time),13,K);
delta_hat_diff_store = zeros(length(time),K);
psi_store = zeros(length(time),13,K);


% 控制器設計與模擬
for k = 1:K
    
    % 初始條件
    x1_0 = 0.25 + 0.02 * rand(1);   % 初始關節角度 (rad)
    x2_0 = 0.33 + 0.02 * rand(1);   % 初始角速度 (rad/s)
    x0 = [x1_0; x2_0]; % 初始狀態
    
   
    x(:, 1) = x0;
    
    T = time;

    xd = [xd_1(k,T);xd_2(k,T);xd_3(k,T)];
    xd_store(:,:,k) = xd(:,:);
   
    e_10 = x(1, 1) - xd(1, 1);
    e_20 = x(2, 1) - xd(2, 1);

    for t = 1:length(T)

        disp(t)
        disp(k)

        % 計算誤差
        % 過渡軌跡公式 (23)
        e_1 = x(1, t) - xd(1, t);
        e_2 = x(2, t) - xd(2, t);
       
        syms tau 
        e_ref_t = @(tau) ((e_10+e_20*tau) * (6 * (1 - tau / Td).^5 - 15 * (1 - tau / Td).^4 + 10 * (1 - tau / Td).^3));
        e1_star_t = piecewise(tau >= Td, 0, tau < Td, e_ref_t);
     
        e1_star = subs(e1_star_t,tau,t*dt);
        e1_star = double(e1_star);
    
        e2_star_t = diff(e1_star_t,tau,1);
        e2_star_diff_t = diff(e1_star_t,tau,2);

        e2_star = subs(e2_star_t,tau,t*dt);
        e2_star = double(e2_star);
        e2_star_diff = subs(e2_star_diff_t,tau,t*dt);
        e2_star_diff = double(e2_star_diff);

        
        e1_store(t,k) = e_1;
        e2_store(t,k) = e_2;
        e1_star_store(t,k) = e1_star;
        e2_star_store(t,k) = e2_star;
        e2_star_diff_store(t,k) = e2_star_diff;
        
        z_1 = e_1 - e1_star;
        z_2 = e_2 - e2_star;
        z = [z_1;z_2];

        s = [c_1,1] * z;
        x_bar = [x(1, t);x(2, t);xd(1, t);xd(2, t);xd(3, t);e1_star;e2_star;e2_star_diff];
        
        rho = 1 / (((kc_upper-z_1)^2) * ((kc_lower+z_1)^2));
        


       
        for j = 1:length(L_centers)

            L = L_centers(j);

            for i = 1:length(x_bar)
                
                X = x_bar(i);

                % 自適應模糊控制器
                omega = psi_fls_M1(X,L); % 激活函數
                
                if i == 1
                    omega_sum = 1;
                end

                omega_sum = omega_sum * omega;
            end
            
            phi(j) = omega_sum;

            if j == 1
                phi_sum = 0;
            end

            phi_sum = phi_sum + phi(j);
        end

        psi(t,:) = phi / phi_sum;
        
        if k == 1
            w_hat(1,:,k) = zeros(1,13); % 初始權重假設
            
            w_hat_diff(t,:) = (-gamma_1/1-gamma_1) * w_hat(t,:,k) + (beta1 * rho * s) * psi(t,:) / (1-gamma_1);

           
                
        else
            w_hat(1,:,k) = w_hat(length(time),:,k-1);
            w_hat_diff(t,:) = (-gamma_1/1-gamma_1) * w_hat(t,:,k) + (gamma_1/1-gamma_1) * w_hat(t,:,k-1) + (beta1 * rho * s) * psi(t,:) / (1-gamma_1);

           

        end

        if k == 1
            delta_hat(1,k) = zeros(1);
            %delta_hat(:,k) = zeros(length(time),1);
            delta_hat_diff(t,k) = (-gamma_2/1-gamma_2) * delta_hat(t,k) + beta2 * (abs(rho*s)) / (1-gamma_2);
            
            

        else
            delta_hat(1,k) = delta_hat(length(time),k-1);
            delta_hat_diff(t,k) = (-gamma_2/1-gamma_2) * delta_hat(t,k) + (gamma_2/1-gamma_2) * delta_hat(t,k-1) + beta2 * (abs(rho*s)) / (1-gamma_2);
            
           

        end

        w_hat(t+1,:,k) = w_hat(t,:,k) + w_hat_diff(t,:) * dt;
        delta_hat(t+1,k) = delta_hat(t,k) + delta_hat_diff(t,k) * dt;
        
        u(t) = -lambda_1 * s / rho - delta_hat(t,k) * sign(rho*s) - w_hat(t,:,k) * psi(t,:)';
        if u(t) >= sp
            u(t) = u_max;
        
        elseif dp <= u(t) && u(t) < sp
            u(t) = (sp - 2 + 2 * sin((5*pi*u(t))/(2*sp))) * ((u(t)-dp)/(sp-dp));
        
        elseif dn <= u(t) && u(t) < dp
            u(t) = 0;

        elseif sn <= u(t) && u(t) < dn
            u(t) = (sn + 2 - 2 * sin((5*pi*u(t))/(2*sn))) * ((u(t)-dn)/(sn-dn));
        
        else
            u(t) = u_min;
        
        end

        dx2 = (m * g * l / J) * sin(x(1, t)) + (1 / J) * u(t);
        x(2, t+1) = x(2, t) + dx2 * dt;
        dx1 = x(2, t);
        x(1, t+1) = x(1, t) + dx1 * dt;
        
         
    end

    psi_store(:,:,k) = psi(:,:);
    w_hat_diff_store(:,:,k) = w_hat_diff(1:length(T),:);
    delta_hat_diff_store(:,k) = delta_hat_diff(t);
    u_store(:,:,k) = u(1,:);
    x1_store(:,:,k) = x(1, 1:length(T));
    x2_store(:,:,k) = x(2, 1:length(T));
end

% 繪製結果
for i = 1:10

figure;
subplot(3, 2, 1);
plot(time, x1_store(:,:,i), 'r', 'LineWidth', 1.5); hold on;
plot(time, xd_store(1,:,i), 'b--', 'LineWidth', 1.5);
xlabel('Time (s)'); ylabel('Joint Angle (rad)');
legend('Actual Angle', 'Reference Angle');
title(['第',num2str(i),'次迭代','Tracking Performance of Joint Angle']);

subplot(3, 2, 2);
plot(time, x2_store(:,:,i), 'r', 'LineWidth', 1.5); hold on;
plot(time, xd_store(2,:,i), 'b--', 'LineWidth', 1.5);
xlabel('Time (s)'); ylabel('Angle speed (rad/s)');
legend('Actual Angle speed', 'Reference Angle speed');
title(['第',num2str(i),'次迭代','Tracking Performance of Angle speed']);

subplot(3, 2, 3);
plot(time, e1_store(:,i), 'r', 'LineWidth', 1.5); hold on;
plot(time, e2_store(:,i), 'b--', 'LineWidth', 1.5); 

subplot(3, 2, 4);
plot(time, w_hat_diff_store(:,5,i));
title('w_hat_diff');

subplot(3, 2, 5);
plot(time, delta_hat_diff_store(:,i));
title('delta_hat_diff');

subplot(3, 2, 6);
plot(time, u_store(:,:,10));
title('u');
end


