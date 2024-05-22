clc;
clear;
tic;
a = 0; %上界
b = pi / 2; %下界
m = 100; %种群规模
d = 2; %变量个数
maxgen = 2000; %最大迭代次数
pc = 0.8; %交叉概率
pm = 0.1; %变异概率
sigma = 0.1; %变异参数
q = 1; %变异点数 这里表示随机两个数变异
elite_ratio = 0.1; %精英保留比例

%初始化种群
pop = repmat(a, m, d) + repmat((b-a), m, d) .* rand(m, d); %生成m*d阶随机矩阵；

history_pop = zeros(m, d, maxgen); % 存储每次迭代的种群分布
history_z = zeros(maxgen, 1); % 存储每次迭代的最大函数值

x = pop(:, 1);
y = pop(:, 2);

for it = 1:maxgen
    
    z = fun(x, y);
    f = z - min(z) + 0.01; % 适应度函数，减去最小值并加上一个小常数，使函数值总是大于零，保持单调性
    
    % 精英保留策略
    elite_num = floor(m * elite_ratio);
    [~, elite_idx] = sort(f, 'descend');
    elite_pop = pop(elite_idx(1:elite_num), :);
    
    %进行轮盘赌
    s = rws1(f);
    pop = pop(s, :); %因为s是值的序号，复制到pop中
    
    %进行交叉
    pop = blendingcrossover(pop , pc);
    
    %进行变异
    pop = mutation_real(pop, sigma, q, pm);
    
    % 更新种群
    pop(1:elite_num, :) = elite_pop;
    
    x = pop(:, 1);
    y = pop(:, 2);
    
    z = fun(x, y);
    [~, v] = sort(z, 'descend'); % 排序顺序修改为降序
    z = z(v);
    history_z(it) = max(z);
    disp(['迭代次数：', num2str(it),'最大函数值为：', num2str(max(z))])
end

figure(1)
plot3(x, y, z, 'o');
xlabel('x');
ylabel('y');
zlabel('z');
title('种群分布图');

figure(2)
plot(history_z);
xlabel('迭代次数');
ylabel('最大函数值');
title('函数值');


function z=fun(x,y)
    z=x/pi+y/pi;%求z的最大值
    0<x<pi/2;0<y<pi/2;
end

function s=rws1(f)
    %轮盘赌选选择；
    %适应度 f;
    %选择结果（个体序号） s;

    n=length(f);
    p=f./sum(f);%适应度值归一化
    ps=cumsum(p);%适应度值累加
    r=rand(n,1);%产生n个随机数
    s=zeros(n,1);
    for i=1:n;
        s(i)=find(r(i)<ps,1);%查找ps(j-1)<r<ps(j) 的j，使得r(i)<ps(j)
    end
end

function pop=blendingcrossover(pop,pc) %混合交叉

    [H,L]=size(pop);
    if (L<3),pc=0.9; end

    k=randperm(H);%生成1~H的随机排列，用于随机配对
    m=100;
    for i=1:m
        if (rand()<pc &&i~=k(i)) %以概率pc进行交叉
            u=rand(1,L);%生成交叉的随机数
            pop(i,:)=u.*pop(i,:)+(1-u).*pop(k(i),:);%混合交叉
        end
    end
end

function pop = mutation_real(pop,sigma, q, pm) %高斯型变异

    %sigma 参数，控制变异程序
    %q 变异点数
    %pm 变异概率
    a = 0; %上界
    b = pi / 2; %下界

    [~,L]=size(pop);
    if (L<4),sigma=0.1;end
    if (L<5),q=1;end
    if (L<6),pm=0.1;end
    m=100;
    for i=1:m
        if rand()<pm %以0.1的概率进行变异
            r=floor(1+L*rand(1,q)); %生成随机变异位数，floor函数，朝负无穷大方向取整
            temp=pop(i,:);
            temp(1,r)=pop(i,r)+sigma*randn(1,q);%randn为正态分布伪随机数 temp函数，把函数中的所以数都拿出来，一 一排列

            if (all(a<=temp & temp<=b)),pop(i,:)=temp;
            end
        end
    end
end
