clear;
clc;
close all;
Clist=[1304 2312;3639 1315;4177 2244;3712 1399;3488 1535;3326 1556;3238 1229;...
    4196 1044;4312  790;4386  570;3007 1970;2562 1756;2788 1491;2381 1676;...
    1332  695;3715 1678;3918 2179;4061 2370;3780 2212;3676 2578;4029 2838;...
    4263 2931;3429 1908;3507 2376;3394 2643;3439 3201;2935 3240;3140 3550;...
    2545 2357;2778 2826;2370 2975];           % 位置坐标 （x, y）

CityNum=size(Clist,1);                        % 问题的规模

dislist=zeros(CityNum); 

for i = 1:CityNum                             % 产生一个大小为 CityNum x CityNum 的 dislist
    for j = 1:CityNum
        dislist(i,j) = ((Clist(i, 1) - Clist(j, 1))^2 + (Clist(i, 2) - Clist(j, 2))^2)^0.5;       
    end
end

TabuList = zeros(CityNum);                      % 禁忌表
TabuLength = round((CityNum * (CityNum - 1) / 2)^0.5);% 禁忌表长度
Candidates = 200;                               % 候选集的个数 (全部领域解个数)
CandidateNum = zeros(Candidates,CityNum);       % 候选解集合
S0 = randperm(CityNum);                         % 随机产生初始解
BSF = S0;                                       % best so far 渴望水平函数：当前最优解
BestL = Inf;                                    % 当前最佳解距离
p=1;                                            % 记录迭代次数
StopL = 500;                                    % 最大迭代次数

figure(1);

if Candidates > CityNum * (CityNum - 1)/2
    disp('候选解个数不大于n*(n-1)/2!');
end
ArrBestL = zeros(1, StopL);
ALong = zeros(1, StopL);
for p = 1 : 1 : StopL

    ALong(p) = Fun(dislist,S0);       % 当前解的适配值
    i = 1;

    A = zeros(Candidates,2);          % 解中交换的城市矩阵

    % 以下while的 是生成随机的200 X 2  的矩阵矩阵A。每一个元素都是在1-31之间的

    while i <= Candidates        
        M = CityNum * rand(1,2);
        M = ceil(M);
        
        if M(1)~=M(2)
            A(i,1)=max(M(1),M(2));        % 选择要进行位置交换的两个城市
            A(i,2)=min(M(1),M(2));
            isSame = 0;                   % 当前产生两个随机值是否之前产生过 1 是 0 否
            if i > 1
                for j = 1:i - 1
                    if A(i,1) == A(j,1) && A(i,2) == A(j,2)
                        isSame = 1;
                        break;
                    end
                end
            end 

            if isSame == 0
               i = i + 1;
            end            
        end
    end

    %---------------- 产生领域解 ----------------------%

    BestCandidateNum = 100;     % 保留前较好的领域解

    BestCandidate = Inf * ones(BestCandidateNum,4);

    F=zeros(1,Candidates);

    % 产生一个S0的邻域
    for i = 1 : Candidates
        CandidateNum(i,:) = S0;  % 候选解集合
        CandidateNum(i,[A(i,2),A(i,1)]) = S0([A(i,1),A(i,2)]);
        F(i) = Fun(dislist,CandidateNum(i,:));

        if i <= BestCandidateNum
            BestCandidate(i,2) = F(i);            % 当前候选解的函数适配值
            BestCandidate(i,1) = i;               % 当前候选解使用的交换规则 A 的序号
            BestCandidate(i,3) = S0(A(i,1));      % 交换的两个城市号
            BestCandidate(i,4) = S0(A(i,2));   
        else
            for j = 1:BestCandidateNum
                if F(i) < BestCandidate(j,2)      % 搜索集的大小为100，而交换集的大小为200，保留其中适配值最优的100个
                    BestCandidate(j,2) = F(i);
                    BestCandidate(j,1) = i;
                    BestCandidate(j,3) = S0(A(i,1));
                    BestCandidate(j,4) = S0(A(i,2));
                    break;
                end            
            end
        end
    end

    

    [JL,Index] = sort(BestCandidate(:,2));    % 对邻域按适应值进行排序
    SBest = BestCandidate(Index,:);
    BestCandidate = SBest;

	if BestCandidate(1,2) < BestL
        BestL = BestCandidate(1,2);            % 当前邻域内适应值优于历史最优值，不考虑禁忌表
        S0 = CandidateNum(BestCandidate(1,1),:);        
        BSF = S0;
        for m = 1 : CityNum
            for n = 1 : CityNum
                if TabuList(m,n) ~= 0
                    TabuList(m,n) = TabuList(m,n) - 1;                  % 更新禁忌表
                end
            end
        end
        TabuList(BestCandidate(1,3),BestCandidate(1,4)) = TabuLength;   % 更新禁忌表
    else  
        for i = 1:BestCandidateNum
            if  TabuList(BestCandidate(i,3),BestCandidate(i,4)) == 0
                S0 = CandidateNum(BestCandidate(i,1),:);                % 禁忌搜索的关键性语句              
                for m=1:CityNum
                    for n=1:CityNum
                        if TabuList(m,n)~=0
                            TabuList(m,n)=TabuList(m,n)-1;               
                        end
                    end
                end        
                TabuList(BestCandidate(i,3),BestCandidate(i,4))=TabuLength; % 更新禁忌表
                break;

            end
        end
	end

    ArrBestL(p)=BestL; 

    for i=1:CityNum-1
        plot([Clist(BSF(i),1),Clist(BSF(i+1),1)],[Clist(BSF(i),2),Clist(BSF(i+1),2)],'bo-');
        hold on;
    end

    plot([Clist(BSF(CityNum),1),Clist(BSF(1),1)],[Clist(BSF(CityNum),2),Clist(BSF(1),2)],'ro-');
    title(['迭代次数:',int2str(p),' 优化最短距离:',num2str(BestL)]);
    hold off;
    pause(0.005);

end

BestShortcut = BSF;                            %最佳路线
theMinDistance = BestL;                        %最佳路线长度

figure(2);
plot(ArrBestL,'b');
xlabel('迭代次数');
ylabel('目标函数值');
title('适应度进化曲线');
grid;
hold on;
plot(ALong,'r')
 
