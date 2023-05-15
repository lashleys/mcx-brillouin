%% part 1
%                     tmp0=(1.f-prop.g*prop.g)/(1.f-prop.g+2.f*prop.g*rand_next_zangle(t));
%                     tmp0*=tmp0;
%                     tmp0=(1.f+prop.g*prop.g-tmp0)/(2.f*prop.g);
% 
%                                         tmp0=fmax(-1.f, fmin(1.f, tmp0));
% 
%                     theta=acosf(tmp0);

%                     #define logistic_uniform(v)  (acos(1.0-2.0*(v))*R_PI)

%                     __device__ RandType rand_uniform01(RandType t[RAND_BUF_LEN])
%                         rand_need_more(t);
%                         return logistic_uniform(t[0]);

% rand_next_zangle(t) = acos(1 - 2*t)/pi
% 
% __device__ void logistic_step(RandType *t, RandType *tnew){
%     t[0]=FUN(t[0]);
%     t[1]=FUN(t[1]);
%     t[2]=FUN(t[2]);
%     t[3]=FUN(t[3]);
%     t[4]=FUN(t[4]);
%     tnew[4]=RING_FUN(t[0],t[4],t[1]);   /* shuffle the results by separation of 1*/
%     tnew[0]=RING_FUN(t[1],t[0],t[2]);
%     tnew[1]=RING_FUN(t[2],t[1],t[3]);
%     tnew[2]=RING_FUN(t[3],t[2],t[4]);
%     tnew[3]=RING_FUN(t[4],t[3],t[0]);
% }
% 
% #define RING_FUN(x,y,z)        (NU2*(x)+NU*((y)+(z)))
% 
% 
%   #define NU                   1e-14
%   #define NU2                  (1.0-2.0*NU)
% 
%     #define FUN(x)               (4.0*(x)*(1.0-(x)))
% 
%       typedef float RandType;
%   #define FUN(x)               (4.f*(x)*(1.f-(x)))
%   #define NU                   1e-7f
%   #define NU2                  (1.f-2.f*NU)
% 
for i = 1:1000000
% __device__ RandType rand_uniform01(RandType t[RAND_BUF_LEN]){
    NU = 1e-7;
    NU2 = 1 - 2*NU;
    g = 0.9;
    % RandType tnew[RAND_BUF_LEN]={0.f};
    t = [rand(1), rand(1), rand(1), rand(1), rand(1)];
    tnew = [0 0 0 0 0];
    t(1)=4*t(1)*(1 - t(1));
    t(2)=4*t(2)*(1 - t(2));
    t(3)=4*t(3)*(1 - t(3));
    t(4)=4*t(4)*(1 - t(4));
    t(5)=4*t(5)*(1 - t(5));
    tnew(5)=NU2*t(1)+NU*(t(5)+t(2));   % /* shuffle the results by separation of 1*/
    tnew(1)=NU2*t(2)+NU*(t(1)+t(3));
    tnew(2)=NU2*t(3)+NU*(t(2)+t(4));
    tnew(3)=NU2*t(4)+NU*(t(3)+t(5));
    tnew(4)=NU2*t(5)+NU*(t(4)+t(1));
    
    tnew(1)=4*tnew(1)*(1 - tnew(1));
    tnew(2)=4*tnew(2)*(1 - tnew(2));
    tnew(3)=4*tnew(3)*(1 - tnew(3));
    tnew(4)=4*tnew(4)*(1 - tnew(4));
    tnew(5)=4*tnew(5)*(1 - tnew(5));
    t(5)=NU2*tnew(1)+NU*(tnew(5)+tnew(2));   % /* shuffle the results by separation of 1*/
    t(1)=NU2*tnew(2)+NU*(tnew(1)+tnew(3));
    t(2)=NU2*tnew(3)+NU*(tnew(2)+tnew(4));
    t(3)=NU2*tnew(4)+NU*(tnew(3)+tnew(5));
    t(4)=NU2*tnew(5)+NU*(tnew(4)+tnew(1));
    rand_next_zangle = acos(1 - 2*t(1))/pi;
                    
    tmp0=(1 - g*g)/(1-g+2*g*rand_next_zangle);
    tmp0 = tmp0*tmp0;                              
    tmp0=(1+g*g-tmp0)/(2*g);                                 
    tmp0=max(-1, min(1, tmp0));
    cos_theta_alt(i) = tmp0;
%     theta(i)=acos(tmp0);
end
figure(1);
histogram(cos_theta_alt,'BinWidth',0.01)
set(gca, 'YScale','log')
axis([-1 1 10 10^6])
%     alpha = 
%     zeta = (1 + g^2 - (1 - g^2 - ((1 - g^2)/(1 - g + 2*g*alpha))^2))/(2*g);
%% part 2
    %close all
    theta_hg = linspace(-pi,pi,1000);
    p_cos_theta = 0.5*(1 - g^2)./(1 + g^2 - 2.*g.*cos(theta_hg)).^(3/2);
    figure(2);
    plot(theta_hg,p_cos_theta)
    set(gca,'YScale','log')
    
    alpha = rand(1,1000000);
%     cos_theta = (1 + g.^2 - (1 - g.^2 - ((1 - g.^2)./(1 - g + 2.*g.*alpha)).^2))./(2.*g);
    cos_theta = ((1 + g.^2 - ((1 - g.^2)./(1 - g + 2.*g.*alpha)).^2))./(2.*g);
%     cos_theta = ()./(2.*g);
%     p_cos_theta = 0.5*(1 - g^2)./(1 + g^2 - 2.*g.*cos_theta).^(3/2);
    figure(3);
    histogram(cos_theta, 'BinWidth',0.01);
%     plot(cos_theta,p_cos_theta)
    set(gca,'YScale','log')
    axis([-1 1 10 10^6])
% RAND_BUF_LEN = 5;
% RandType t[RAND_BUF_LEN];
% gpu_rng_init(t,n_seed,idx);
% 
% __device__ void gpu_rng_init(RandType t[RAND_BUF_LEN],uint *n_seed,int idx){
%     logistic_init(t,n_seed,idx);
% 
% __device__ void logistic_init(RandType *t,uint seed[],uint idx){
%      int i;
%      for(i=0;i<RAND_BUF_LEN;i++)
%            t[i]=(RandType)seed[idx*RAND_BUF_LEN+i]*R_MAX_C_RAND;
%  #define R_MAX_C_RAND         (1.f/RAND_MAX)
%  idx= blockDim.x * blockIdx.x + threadIdx.x;
%  // rand_max >= 32767
%      for(i=0;i<INIT_LOGISTIC;i++)  /*initial randomization*/
%            rand_need_more(t);
% }
% 
% __device__ void rand_need_more(RandType t[RAND_BUF_LEN]){
%     RandType tnew[RAND_BUF_LEN]={0.f};
%     logistic_step(t,tnew);
%     logistic_step(tnew,t);
% }
% 
% __device__ void logistic_step(RandType *t, RandType *tnew){
%     t[0]=FUN(t[0]);
%     t[1]=FUN(t[1]);
%     t[2]=FUN(t[2]);
%     t[3]=FUN(t[3]);
%     t[4]=FUN(t[4]);
%     tnew[4]=RING_FUN(t[0],t[4],t[1]);   /* shuffle the results by separation of 1*/
%     tnew[0]=RING_FUN(t[1],t[0],t[2]);
%     tnew[1]=RING_FUN(t[2],t[1],t[3]);
%     tnew[2]=RING_FUN(t[3],t[2],t[4]);
%     tnew[3]=RING_FUN(t[4],t[3],t[0]);
% }