%% ANÁLISE DINÂMICA 
%
% TIPO DO MODELO DA ESTRUTURA
% Modelo de Shear Building base fixa : tm = 1;
% Modelo de Shear Building com ISE : tm = 2; 
% Shear Building equipado com PTLCD: tm = 3;
% Shear Building equipado com PTLCD e ISE: tm = 4;
  tm = 4;
% TIPO DE SOLO
% Solo denso: s = 1;
% Solo médio: s = 2;
% Solo fofo:  s = 3;
  s = 3;
% TIPO DE ANÁLISE
% Frequencia vs Resposta estrutural: ta = 1 ;
% Histórico da resposta dinâmica: ta = 2 ;
  ta = 2;   
% TIPO DE CARREGAMENTO:
% Aplicado aos nós da estrutura: tc = 1;  
% Aceleração na base - força efetiva nos nós da estrutura : tc = 2;
% Sísmo gerado a partir do histórico de acelerações: tc = 3;
  tc = 3;
% NUMERAÇÃO SISMOS
% El Centro - sismo 1   (dt = 0,02 s)
% Northridge - sismo 2  (dt = 0,005 s)
% Kobe - sismo 3        (dt = 0,01 s)
  sismo = 3; 
% TIPO DE GEOMETRIA DA FUNDAÇÃO
% Circular: f = 1;
% Retangular: f = 2;
  gf = 1;
% OPÇÃO DE PLOTAGEM DA RESPOSTA DINÂMICA
% Plotar de resposta em termos de deslocamento: tp = 1;
% Plotar de resposta em termos de velocidade: tp = 2;
% Plotar de resposta em termos de aceleração: tp = 3;
  tp = 1;
% 
%% PARÂMETROS DA ESTRUTURA 
%  
    format long;
%  
% PROPRIEDADES DA ESTRUTURA
    nel = 40 ;                       % Quantidade de pavimentos
    nnel = 2 ;                       % Quantidade de nós por pilar
    ndof = 1 ;                       % Quantidade de graus de liberdade por nó
    nnode = (nnel-1)*nel+1 ;         % Quantidade total de nós 
    sdof = nnode*ndof ;              % Graus de liberdade do sistema
%     p = 1 ;                         % Quantidade de pilares por pavimento
%     el = 28e+09 ;                   % Modulo de elasticade do material do pilar
%     area = 0.3*0.3 ;                % Seçao transversal do pilar quadrado
%     xi = 6.75e-04 ;                 % Momento de inercia da seçao transversal
    drm = 0.0343 ;                   % Damping ratio m 
    drn = 0.0343 ;                   % Damping ratio n    
    Mi = 9.8E+5 ;                    % Massa concentrada no pavimento
    Iye = 1.31E+08 ;                 % Momento massa de inercia do pavimento
    tleng = 4*nel ;                  % Comprimento total do Shear Building
    leng = tleng/nel ;               % Uniforme mesh(elementos iguais)
%     nbc= 1 ;                        % Quantidade de graus com restriçoes
%     bcdof(1) = 1 ;                  % Grau de liberdade com restrição
%     kk = zeros(sdof,sdof);          % Initialization of system stiffness matrix
%     index = zeros(nnel*ndof,1);     % Initialization of index vector
%
%% PARAMETROS SOLO/FUNDAÇAO
% Solo
%     if s == 1
%      'Stiff Soil : S-1' %#ok<NOPTS>
%         Gs = 1.71E+08;  % Modulo de elasticidade transversal
%         vs = 0.48;      % Coeficiente de Poisson
%         Vs = 300.0;     % Velocidade da onda de cisalhamento
%         rhos = 1900;    % Massa específica do solo 
%         
%      elseif s == 2      
%      'Soft Soil : S-2'  %#ok<NOPTS>
%         Gs = 1.8E+07;   % Modulo de elasticidade transversal 
%         vs = 0.49;      % Coeficiente de Poisson
%         Vs = 100.0;     % Velocidade da onda de cisalhamento 
%         rhos = 1800;    % Massa  específica do solo
%         
%      elseif s == 3      
% 'Very Soft Soil : S-3'  %#ok<NOPTS> 
%         Gs = 5.4E+06;   % Modulo de elasticidade transversal 
%         vs = 0.4;       % Coeficiente de Poisson
%         Vs = 57.3;      % Velocidade da onda de cisalhamento 
%         rhos = 1600;    % Massa específica do solo
%         
%     end
%    
% Fundação    
    if gf == 1                       % Fundaçao circular                                
%         r = 10.61 ;                  % Raio da fundação
%         e = 1.0;                     % Espessura
%         rho = 2500;                  % Massa específica da fundação
%         Mo = (pi*r^2)*e*rho;         % Massa da fundaçao 
%         Io = Mo*(pi*r^2)/4;          % Momento massa da fundaçao 
%         Ks = (8*Gs*r)/(2-vs);           % Rigidez translaçao (x)
%         Kr = (8*Gs*r^3)/(3*(1-vs));     % Rigidez rotaçao (yy)
%         Cs = 4.6*(r^2)*rhos*Vs/(2-vs);  % Amortecimento solo translaçao
%         Cr = 0.4*(r^4)*rhos*Vs/(1-vs);  % Amortecimento solo rotaçao
%
    % Exemplo de validação (Farshidianfar e Soheili (2013))       
        if s == 1
           %'Dense Soil : S-1'  %#ok<NOPTS>
            Ks = 5.75E+10;  % Rigidez translaçao (x)
            Kr = 1.91E+13;  % Rigidez rotaçao (yy)
            Cs = 1.32E+09;  % Amortecimento solo translaçao
            Cr = 1.15E+11;  % Amortecimento solo rotaçao 
         elseif s == 2
           %'Medium Soil : S-2' %#ok<NOPTS>
            Ks = 1.80E+10;  % Rigidez translaçao (x)
            Kr = 7.02E+12;  % Rigidez rotaçao (yy)
            Cs = 6.90E+08;  % Amortecimento solo translaçao
            Cr = 7.02E+10;  % Amortecimento solo rotaçao
         
        elseif s == 3
           %'Soft Soil : S-3'   %#ok<NOPTS>
            Ks = 1.91E+09;  % Rigidez translaçao (x)
            Kr = 7.53E+11;  % Rigidez rotaçao (yy)
            Cs = 2.19E+08;  % Amortecimento solo translaçao
            Cr = 2.26E+10;  % Amortecimento solo rotaçao
            
        end
        Mo = 1.96E+06;  % Massa da fundaçao 
        Io = 1.96E+08;  % Momento massa da fundaçao 
        
     elseif gf == 2                                                                              % Fundação retangular   
%         l = 3;                                                                                 % Metade da dimensão em x  
%         b = 3;                                                                                 % Metade da dimensão em y 
%         R = l/b;                                                                                 % Relação entre as dimensões
%         e = 0.5;                                                                                   % Espessura
%         rho = 2500;                                                                              % Massa específica da fundação
%         Mo = (2*l*2*b)*e*rho;                                                                    % Massa da fundaçao
%         Io = (Mo*(2*l)^2)/6;                                                                     % Momento massa da fundaçao 
%         Ks = ((2*Gs*l)/(2-vs))*(2+2.5*(1/R)^0.85);                                               % Rigidez translaçao
%         Kr = ((Gs)/(1-vs))*(Io^0.75)*((R)^0.25)*(2.4+0.5*(1/R));                                 % Rigidez rotaçao
%         ao = (2*pi*fse(1)*b)/Vs ;                                                               % Frequencia adimensional
%         alphas = 1;                                                                              % Modificadores de rigidez dinâmica
%         alphar = 1-((0.55+(0.01*sqrt(R-1)))*ao^2)/((2.4-(0.4/R^3))+ao^2);  
%         phi = sqrt(2*(1-vs)/(1-2*vs));
%         Bs = (4*R/(Ks/(Gs*b)))*(ao/2);                                                           % Radiation Damping Ratios
%         Br = (((4*phi/3)*(R^3)*(ao^2))/((Kr/(Gs*b^3))*((2.2-0.4/(R^3))+ao^2)))*(ao/(2*alphar));
%         Cs = 2*Bs*Ks/(2*pi*fse(1));                                                             % Amortecimento solo translaçao
%         Cr = 2*Br*Kr/(2*pi*fse(1));                                                             % Amortecimento solo rotaçao
    end    
%
%
%% PARAMETROS DO PTLCD
    % Parametros do TLCD
    L = 40;                       % Comprimento total do tubo (m)
    B = 0.9*L;                    % Comprimento horizontal do tubo (m)  
    D = 1.769;                    % Diametro do tubo (m)
    Aa = pi*((D^2)/4);            % Area do tubo do ACL  
    rhof = 998.207;               % Massa especifica do fluido (kg/m³)
    g = 9.807;                    % Gravidade (m/s²)
    rug = 1.5e-06;                % Rugosidade do tubo (m)
    visc = 1.003e-06;             % Viscosidade cinemática do fluido
    Ra = 0.12;                    % Razão de abertura do orifício (Ao/Aa)
    ck = pi*L*D*rhof/8;                    % Parcela constante do amortecimento 
    cdk = rhof*(Aa/2)*(((1/Ra)-1)^2);      % Parcela de amortecimento causada pelo orifício
    
    Mtlcd = rhof*Aa*L;            % Massa do PTLCD
    Ctlcd = 0;                    % Amortecimento do PTLCD
    
    wtlcd = sqrt((((pi*(D^2)*rhof*g)/2)/Mtlcd));
    ftlcd = wtlcd /(2*pi);
    
    % PARÂMETROS DO TRECHO PRESSURIZADO
    
    % Quantidade de PTLCD sintonizado em dada frequência
    % colocar em forma de vetor coluna, pois cada linha representa um
    % ou mais PTLCD sintonizado em diferentes frequências
        PTLCD_f = [14];
    
    % Quantidade total de atenuadores
        Ntlcd = 0;
        for i = 1:length(PTLCD_f)
            Ntlcd = Ntlcd + PTLCD_f(i); 
        end
    
    % Pressão do gás do PTLCD sintonizado em dada frequência (atm)
    % colocar em forma de vetor coluna, pois cada linha representa uma
    % frequência de funcionamento de cada grupo de PTLCD
        Po = [0.193];  
        
    % Altura do pressurizador (m)
        Z = 2.00;     
    
    Ktlcd = zeros(length(PTLCD_f),1);
    wptlcd = zeros(length(PTLCD_f),1);
    fptlcd = zeros(length(PTLCD_f),1);
    
    for i = 1:length(PTLCD_f)
        % Rigidez do PTLCD   
        Ktlcd(i) = ((pi*(D^2)*rhof*g)/2)+2*1.4*(Po(i)*101325/Z)*(pi*(D^2)/4);
        
        % Frequencia do PTLCD
        wptlcd(i) = sqrt(Ktlcd(i)/Mtlcd);
        fptlcd(i)= wptlcd(i) /(2*pi);  
    end   
%     
%% LOOP PARA A QUANTIDADE TOTAL DE ELEMENTOS
% 
%     for iel=1:nel 
%         index = feeldof1(iel,nnel,ndof);    % Extrair os GDL do sistema para cada elemento
%         [k] = shearbuilding(el,xi,leng,p);  % Computar as matrizes dos elementos
%         kk = feasmbl1(kk,k,index);          % Superposição da matriz de rigidez da estrutura
%     end
%
%% MATRIZ DE RIGIDEZ E MATRIZ DE MASSA CONCENTRADA DO SHEAR BUILDING (Ke, Me)
%
% Matriz de Rigidez do Shear Building
%     gl = length(kk);        % Graus de liberdade 
%     Ke = zeros(gl-1,gl-1);  % Dimensao da matriz de rigidez da estrutura
%     gle = length(Ke);       % Graus de liberdade do Shear Building  
%     for i=2:sdof            % Coeficientes de rigidez do Shear Building (kk)
%        for j=2:sdof
%            Ke(i-1,j-1)=kk(i,j);     
%        end
%     end
%   
    ki = zeros(40,1);
    ki(1)= 2.13E+09;  
    for i = 2:40
        ki(i)= ki(i-1)-2.8769231E+7; %2.8769231e+7
    end
    gle = 40;
    Ke = zeros(gle,gle);
    for i = 1:39
        Ke(i,i) = ki(i,1)+ki(i+1,1);
        Ke(i,i+1) = -ki(i+1,1);
        Ke(i+1,i) = -ki(i+1,1);
    end
    Ke(gle,gle)=ki(gle,1);
% Matriz de Massa Concentrada do Shear Building
    Me = zeros(gle,gle);    % Dimensao da matriz de massa da estrutura
    for i=1:gle             % Coeficientes de massa da estrutura (Mi)
        Me(i,i)= Mi;                    
    end
    Met=(ones(gle,1))'*Me*(ones(gle,1)); % Massa total da estrutura 
%
%% FREQUENCIA, PERIODO e MODOS DE VIBRAÇÃO
    [AVE,AVA]=eig(Ke,Me);           % Problema de autovetor (AVE) e autovalor (AVA)
    num = 1:1:(gle);                % Contador
    wn = zeros(gle,1);              % Inicializaçao do vetor de frequencias naturais (rad/s)
    for i=1:(gle)
        wn(i) = (sqrt(AVA(i,i)));   % Frequencias naturais da estrutura (rad/s)
    end
    fn = zeros(gle,1);              % Inicializaçao do vetor de frequencias naturais (Hz)
    for i=1:(gle)
        fn(i) = wn(i)/(2*pi);       % Frequencias naturais da estrutura (Hz)
    end
    Tn = zeros(gle,1);              % Inicializaçao do vetor de periodos naturais (s)
    for i=1:(gle)
        Tn(i) = (2*pi)/wn(i);       % Periodos naturais da estrutura (s)
    end
    modos = [num' AVE];             % Matriz espectral da estrutura
%  
%% MATRIZ DE AMORTECIMENTO ESTRUTURAL (Ce) 
    Ce = zeros(gle,gle);                         % Matriz de amortecimento       
    [Ce] = propdamping(gle,Ke,Ce,Me,wn,drm,drn); % Matriz de amortecimento proporcional da estrutura 
%
%% PARAMETROS DA ESTRUTURA NO MODO FUNDAMENTAL

    Li = zeros(gle,1);  % Participação modal
    Li2 = zeros(gle,1); % Massa Modal 
    l=0;
    for i=1:gle
        Li(i)= AVE(:,i)'*(Me*(ones(gle,1)));
        Li2(i)= Li(i)^2;
        start = l + Li2(i);
        l = start;
    end
    MME = (1/l)*Li2;                      % Participação modal da massa
    Mef = Met*MME(1);                     % Massa efetiva para o primeiro modo de vibraçao
%     Kest = (4*(pi^2)*Mef)/(tn(1)^2);    % Rigidez lateral do oscilador
%     Cest = 2*drm*sqrt(Kest*Mef);        % Amortecimento do oscilador 
%     Hef = tleng*0.7;                    % Altura efetiva do primeiro modo
%     rk = Hef*wn(1)/Vs;                  % Razao de rigidez entre estrutura e solo
%    
%% INTERAÇÃO SOLO-ESTRUTURA (ISE)
%         
% MATRIZ DE RIGIDEZ ISE (Kse)
    % Coeficientes de rigidez da estrutura (Ke)
        Kse = zeros(gle + 2,gle + 2);   % Dimensao matriz de rigidez ISE
        glse = length(Kse);             % Graus de liberdade do sistema ISE (glse)
        for i=1:gle                     % Rigidez da estrutura
           for j=1:gle
               Kse(i+2,j+2)= Ke(i,j); 
           end
        end
    % Coeficientes de rigidez do solo (Ks)
        Kse(2,2)= Ks;   % Rigidez solo translaçao
        Kse(1,1)= Kr;   % Rigidez solo rotaçao
%      
% MATRIZ DE MASSA ISE (Mse)
    % Matriz de massa estrutura (Me)
        Mse=zeros(glse,glse);   % Dimensao matriz de massa ISE
        for i=1:gle             % Coeficientes de massa da estrutura (Mi)
           for j=1:gle
               Mse(i+2,j+2)= Me(i,j);
           end
        end
    % Matriz de massa solo-estrutura (Mse,Mes)
        % Coeficiente: Mi
            for i=1:gle   % Massa da estrutura  
                Mse(i+2,2)=Me(i,i); 
                Mse(2,i+2)=Me(i,i);
            end
        % Coeficiente: MiZi
            for i=1:gle  
                Mse(i+2,1)=Me(i,i)*i*leng;
                Mse(1,i+2)=Me(i,i)*i*leng;
            end
    % Matriz de massa do solo (Ms)    
        % Coeficiente: M0 + somatorio(Mi)
            SMx = 0;
            for i=1:gle
                M = Me(i,i);
                SMx = SMx + M;
            end
            Mse(2,2) = Mo + SMx;
        %
        % Coeficiente: somatorio(MiZi)
            SMxH = 0;
            for i=1:gle
                MZ = Me(i,i)*i*leng;
                SMxH = SMxH + MZ;
            end
            Mse(2,1) = SMxH;
            Mse(1,2) = SMxH;
        %
        % Coeficiente: I0 + somatorio(Ii + MiZi²)
            SMyy = 0;
            for i=1:gle
                MI = Me(i,i)*((i*leng)^2)+ Iye;
                SMyy = SMyy + MI;
            end
            Mse(1,1) = Io  + SMyy;
        %
% Sistema em Vibração Livre e sem amortecimento
      wse = sqrt(eig(Kse,Mse));
      fse = (1/(2*pi))*sqrt(eig(Kse,Mse));
      Tse = zeros(glse,1);
      for i = 1 : glse
          Tse(i,1) = 1/fse(i);
      end
%      
% MATRIZ DE AMORTECIMENTO ISE (Cse) 
        Cse = zeros(glse,glse);                % Dimensao matriz de amortecimento ISE
        for i=1:gle                          % Coeficientes de amortecimento estrutural
           for j=1:gle
              Cse(i+2,j+2)= Ce(i,j);   
           end
        end
   % Coeficientes de amortecimento do solo (Cs)
        if gf == 1
            Cse(2,2)= Cs;   % Amortecimento solo translaçao
            Cse(1,1)= Cr;   % Amortecimento solo rotaçao
        else
            ao = (2*pi*fse(1)*b)/Vs ;                                                                % Frequencia adimensional
            alphas = 1;                                                                              % Modificadores de rigidez dinâmica
            alphar = 1-((0.55+(0.01*sqrt(R-1)))*ao^2)/((2.4-(0.4/R^3))+ao^2);  
            phi = sqrt(2*(1-vs)/(1-2*vs));
            Bs = (4*R/(Ks/(Gs*B)))*(ao/2);                                                           % Radiation Damping Ratios
            Br = (((4*phi/3)*(R^3)*(ao^2))/((Kr/(Gs*B^3))*((2.2-0.4/(R^3))+ao^2)))*(ao/(2*alphar));
            Cse(2,2) = 2*Bs*Ks/(2*pi*fse(1));                                                        % Amortecimento solo translaçao
            Cse(1,1) = 2*Br*Kr/(2*pi*fse(1));                                                        % Amortecimento solo rotaçao          
        end 
%
%% CONTROLE DE VIBRAÇÃO (PTLCD)
%         
% MATRIZ DE RIGIDEZ  
    % Rigidez da estrutura
        glf = gle + Ntlcd;      % Graus de liberdade (E + MPTLCD)
        Kf = zeros(glf,glf);    % Dimensao matriz de rigidez 
        for i=1:gle             % Rigidez da estrutura
           for j=1:gle
               Kf(i,j)= Ke(i,j); 
           end
        end
    % Atenuador de Coluna Líquida 
        cont = 0;
        for i = 1:length(PTLCD_f)
            for j = 1:PTLCD_f(i)  
                cont = cont + 1;
                Kf(gle+cont,gle+cont)= Ktlcd(i);                           
            end
        end
%      
% MATRIZ DE MASSA DO SISTEMA 
    % Massa da estrutura
        Mf = zeros(glf,glf);    % Dimensao matriz de massa (E+MPTLCD)
        for i=1:gle-1           % Coeficientes de massa da estrutura (Mi)
           for j=1:gle-1
               Mf(i,j)= Me(i,j);
           end
        end
        Mf(gle,gle)= Mi + Ntlcd*Mtlcd;
    % Atenuador de Coluna Líquida 
        for i = 1:Ntlcd            
            Mf(gle, gle+i) = (B/L)*Mtlcd;
            Mf(gle+i, gle) = (B/L)*Mtlcd;
            Mf(gle+i, gle+i) = Mtlcd;
        end
%      
% MATRIZ DE AMORTECIMENTO DO SISTEMA (Ce) 
        Cf = zeros(glf,glf);    % Dimensao matriz de amortecimento(E+MPTLCD)
     % Amortecimento estrutural
        for i=1:gle             % Coeficientes de amortecimento estrutural
           for j=1:gle
              Cf(i,j)= Ce(i,j);   
           end
        end
     % Atenuador de Coluna Líquida 
        for i = 1:Ntlcd
            Cf(gle+i,gle+i)= Ctlcd;     % Rigidez solo translaçao
        end
%
%% INTERAÇÃO SOLO-ESTRUTURA COM PTLCD
%
% MATRIZ DE RIGIDEZ 
        Kia = zeros(glse+Ntlcd,glse+Ntlcd);   % Dimensao matriz de rigidez (ISE+E+TLCD)
    % Coeficientes de rigidez da estrutura (Kia)        
        glia = length(Kia); % Graus de liberdade do sistema 
        for i=1:gle         % GL da estrutura
           for j=1:gle
               Kia(i+2,j+2)= Ke(i,j); 
           end
        end
    % Coeficientes de rigidez do solo     
        Kia(1,1)= Kr;   % Rotaçao em Y
        Kia(2,2)= Ks;   % Translaçao em X
    % Coeficientes de rigidez do TLCD 
        cont = 0;
        for i = 1:length(PTLCD_f)
            for j = 1:PTLCD_f(i)  
                cont = cont + 1;
                Kia(glse+cont,glse+cont)= Ktlcd(i);                           
            end
        end
%
% MATRIZ DE MASSA  
        Mia=zeros(glia,glia);     % Dimensao matriz de massa (ISE+E+TLCD)    
    %
    % Coeficientes de massa da ISE 
        for i=1:gle-1             
           for j=1:gle-1
               Mia(i+2,j+2)= Me(i,j);
           end
        end
        Mia(glse,glse)= Mi + Ntlcd*Mtlcd;
        
    % Coeficientes de massa do TLCD
        for i = 1:Ntlcd 
            Mia(glse+i,glse) = (B/L)*Mtlcd;
            Mia(glse,glse+i) = (B/L)*Mtlcd;
            Mia(glse+i,glse+i) = Mtlcd;
        end
    %    
    % Matriz de massa solo-estrutura (Mse,Mes)
        % Coeficiente: Mi
            for i=1:gle     % Massa da estrutura  
                Mia(i+2,2)=Me(i,i); 
                Mia(2,i+2)=Me(i,i);
                if i == gle
                    Mia(i+2,2)=Me(i,i)+Ntlcd*Mtlcd; 
                    Mia(2,i+2)=Me(i,i)+Ntlcd*Mtlcd;
                end
            end
            for i = 1:Ntlcd
                Mia(glse+i,2)=(B/L)*Mtlcd;
                Mia(2,glse+i)=(B/L)*Mtlcd;
            end
        %    
        % Coeficiente: Mihi
            for i=1:gle  
                Mia(i+2,1)=Me(i,i)*i*leng;
                Mia(1,i+2)=Me(i,i)*i*leng;
                if i == gle
                    Mia(i+2,1)=(Me(i,i)+Ntlcd*Mtlcd)*i*leng; 
                    Mia(1,i+2)=(Me(i,i)+Ntlcd*Mtlcd)*i*leng;
                end
            end
            for i=1:Ntlcd
                Mia(glse+i,1)=(B/L)*Mtlcd*tleng;
                Mia(1,glse+i)=(B/L)*Mtlcd*tleng;
            end
        %    
    % Matriz de massa do solo (Ms)    
        % Coeficiente: M0 + somatorio(Mi)
            SMx = 0;
            for i=1:gle
                M = Me(i,i) + SMx;
                SMx = M;
            end
            Mia(2,2) = Mo + SMx + Mtlcd;
        %
        % Coeficiente: somatorio(MiZi)
            SMxH = 0;
            for i=1:gle
                MZ = Me(i,i)*i*leng;
                SMxH = SMxH + MZ;
            end
            Mia(2,1) = SMxH + Ntlcd*Mtlcd*tleng;
            Mia(1,2) = SMxH + Ntlcd*Mtlcd*tleng;
        %    
        % Coeficiente: I0 + somatorio(Ii + MiZi²)
            SMyy = 0;
            for i=1:gle
                MI = Me(i,i)*((i*leng)^2) + Iye;
                SMyy = SMyy + MI;
            end
            Mia(1,1) = Io + SMyy + Ntlcd*Mtlcd*(tleng)^2;
        % 
% MATRIZ DE AMORTECIMENTO DO SISTEMA (Cia)
    Cia = zeros(glia,glia); % Dimensao matriz de amortecimento
    % Coeficientes de amortecimento estrutural
        for i=1:gle             
           for j=1:gle
              Cia(i+2,j+2)= Ce(i,j);   
           end
        end
     % Coeficientes de amortecimento do solo  
        Cia(1,1) = Cse(1,1);    % Rotação
        Cia(2,2) = Cse(2,2);    % Translação
    % Coeficientes de rigidez do TLCD 
        Cia(glia,glia) = Ctlcd;
%
%% RESPOSTA NO DOMINIO DA FREQUENCIA
%
if ta == 1
        dt = 0.01;              % Time step size
        ti = 0;                 % Initial Time
        tf = 60;                % Final time    
        nt = ((tf-ti)/dt);      % Number of time steps
        Po = 1;                 % Amplitude da força externa
        passo = 0.01;           % Passo frequência
        frf = 0.6;              % Frequência final
        intt = frf/passo;       % Intervalo Total
        faf = zeros(intt+1,1);  % Faixa de frequências
        xm = zeros(intt+1,gle); % Vetor com os valores de deslocamento máximo
        vm = zeros(intt+1,gle); % Vetor com os valores de velocidade máxima
        am = zeros(intt+1,gle); % Vetor com os valores de aceleração máxima
    %    
    % Parametros constantes para a integraçao numerica pelo Metodo de Newmark
        a1 = 1/(0.25*(dt^2));
        a2 = 0.5/(0.25*dt);
        a3 = 1/(0.25*dt);
        a4 = 0.5/0.25;
        a5 = 1/(2*0.25);
        a6 = (dt)*((0.5/(2*0.25))-1);
    %
    if tm == 1
      %
      % INICIALIZAÇAO DOS VETORES
        for i = 0:1:intt       
            faf(i+1) = i*passo;     % Vetor com as frequencias de excitaçao (Hz)
            we = (2*pi*i*passo);    % Frequencia de excitaçao (rad/s)
            x = zeros(gle,nt+1);    % Inicializaçao da matriz de deslocamento
            v = zeros(gle,nt+1);    % Inicializaçao da matriz de velocidade
            a = zeros(gle,nt+1);    % Inicializaçao da matriz de aceleraçao
            dx = zeros(gle,1);      % Incremento da matriz de deslocamento
            dv = zeros(gle,1);      % Incremento da matriz de velocidade
            da = zeros(gle,1);      % Incremento da matriz de aceleração
            v(:,1) = zeros(gle,1);  % Velocidade inicial: {v}o = 0
            x(:,1) = zeros(gle,1);  % Deslocamento inicial: {x}o = 0
            P = zeros(gle,nt);      % Inicializaçao da vetor força    
            dFef = zeros(gle,1);    % Incremento do vetor força
            Po = 1;                 % Amplitude do carregamento 
        %
        % CARREGAMENTO APLICADO  
            if tc == 1 
                for j = 0:1:(nt+1)
                    P(:,j+1) = Po*sin(we*(j)*dt);
                end   

            elseif tc == 2
                U = ones(gle,1);
                for j = 0:1:(nt+1)    
                    P(:,j+1) = -Me*U*(sin(we*(j)*dt));
                end
            end
        % Metodo de Newmark
            % Entradas para iniciar o algoritmo
                a(:,1) = Me\(P(:,1) - Ce*v(:,1) - Ke*x(:,1));
                Kef =  a1*Me + a2*Ce+ Ke;
                for k = 1:1:nt
                        dFef = (P(:,k+1)-P(:,k)) + (a3*Me+ a4*Ce)*v(:,k) + (a5*Me+ a6*Ce)*a(:,k) ;
                        dx = Kef\dFef;                    
                        dv = a2*dx - a4*v(:,k) - a6*a(:,k);
                        da = a1*dx - a3*v(:,k)- a5*a(:,k);
                        x(:,k+1) = x(:,k) + dx;
                        v(:,k+1) = v(:,k) + dv;
                        a(:,k+1) = a(:,k) + da; 
                end
            % Deslocamento máximo
                for l = 1:gle
                    if abs(max(x(l,:))) > abs(min(x(l,:)))
                        xm(i+1,l) = max(x(l,:));
                    else
                        xm(i+1,l) = abs(min(x(l,:)));
                    end
                end
            % Velocidade máxima
                for l = 1:gle
                    if abs(max(v(l,:))) > abs(min(v(l,:)))
                        vm(i+1,l) = max(v(l,:));
                    else
                        vm(i+1,l) = abs(min(v(l,:)));
                    end                
                end
            % Aceleração máxima
                for l = 1:gle
                    if abs(max(a(l,:))) > abs(min(a(l,:)))
                        am(i+1,l) = max(a(l,:));
                    else
                        am(i+1,l) = abs(min(a(l,:)));
                    end                  
                end
        end
    %     
    elseif tm == 2                    % Modelo shear building com ISE                   
        for i = 0:1:intt              % Passo frequencia 
            faf(i+1) = i*passo;       % Vetor com as frequencias de excitaçao (Hz);        
            we = (2*pi*i*passo);      % Frequencia de excitaçao (rad/s)
            a = zeros(glse,nt+1);     % Inicializaçao da matriz de aceleraçao
            da = zeros(glse,1);       % Incremento da matriz de aceleração
            at = zeros(glse-2,nt+1);  % Inicializaçao da matriz de aceleraçao total
            v = zeros(glse,nt+1);     % Inicializaçao da matriz de velocidade
            dv = zeros(glse,1);       % Incremento da matriz de velocidade
            vt = zeros(glse-2,nt+1);  % Inicializaçao da matriz de velocidade total
            x = zeros(glse,nt+1);     % Inicializaçao da matriz de deslocamento
            dx = zeros(glse,1);       % Incremento da matriz de deslocamento
            u = zeros(glse-2,nt+1);   % Inicializaçao da matriz de deslocamento hprizontal total
            v(:,1) = zeros(glse,1);   % Velocidade inicial: {v}o = 0
            x(:,1) = zeros(glse,1);   % Deslocamento inicial: {x}o = 0
            dFef = zeros(glse,nt);    % Vetor de força efetiva
            P = zeros(glse,nt);       % Inicializaçao da vetor força
            Po = 1;                   % Amplitude do carregamento           
            Ft = zeros(1,nt);         % Inicializaçao do vetor força GL de translaçao do solo  
            Mt = zeros(1,nt);         % Inicializaçao do vetor força GL de rotaçao do solo 
            Meff = zeros(glse,1); 
            U = zeros(glse,1);
            U(2,1) = 1;
        %
        % CARREGAMENTO APLICADO
                if tc == 1            % Força aplicada na estrutura
                    for j = 1:gle               
                      for it = 0:1:nt
                          P(j+2,it+1) = Po*sin(we*(it)*dt);
                          P(1,it+1) = Mt(1,it+1) + P(j+2,it+1)*j*leng;
                          P(2,it+1) = Ft(1,it+1) + P(j+2,it+1);
                      end
                      Mt = P(1,:);  
                      Ft = P(2,:);     
                    end
                %
                elseif tc == 2       % Força aplicada na base da estrutura                    
                    Meff = Mse*U;
                    for it=0:1:nt    
                        P(:,it+1) = -Meff*(sin(we*(it)*dt));
                    end
                end
           % Metodo de Newmark
                % Entradas para iniciar o algoritmo           
                a(:,1) = Mse\(P(:,1) - Cse*v(:,1) - Kse*x(:,1));
                Kef =  a1*Mse + a2*Cse+ Kse;
                for it = 1:1:(nt-1)
                        dFef = (P(:,it+1)-P(:,it)) + (a3*Mse+ a4*Cse)*v(:,it) + (a5*Mse+ a6*Cse)*a(:,it) ;
                        dx = Kef\dFef;                    
                        dv = a2*dx - a4*v(:,it) - a6*a(:,it);
                        da = a1*dx - a3*v(:,it) - a5*a(:,it);
                        x(:,it+1) = x(:,it) + dx;
                        v(:,it+1) = v(:,it) + dv;
                        a(:,it+1) = a(:,it) + da; 
                end
            % Resposta do deslocamento horizontal total do sistema acoplado
                %for k = 1:gle
                %    u(k,:)= x(k+2,:) + x(2,:) + k*leng*x(1,:);
                %end
                % Deslocamento máximo {u}
                    for l = 1:gle
                        if abs(max(x(l,:))) > abs(min(x(l,:)))
                            xm(i+1,l) = max(x(l,:));
                        else
                            xm(i+1,l) = abs(min(x(l,:)));
                        end
                    end
            % Resposta absoluta em termos de velocidade do sistema acoplado
                %for k = 1:gle
                %    vt(k,:)= v(k+2,:) + v(2,:) + k*leng*v(1,:);
                %end
                % Velocidade absoluta máxima {vt}
                %    for l = 1:gle
                %        vm(i+1,l)= max(vt(l,:));
                %    end
            % Resposta absoluta em termos de velocidade do sistema acoplado
                %for k = 1:gle
                %    at(k,:)= a(k+2,:) + a(2,:) + k*leng*a(1,:);
                %end
                % Aceleração absoluta máxima {at}
                %    for l = 1:gle
                %        am(i+1,l)= max(at(l,:));
                %    end
        end
    elseif tm == 3
        for i = 0:1:intt              % Passo frequencia 
                faf(i+1) = i*passo;     % Vetor com as frequencias de excitaçao (Hz);        
                we = (2*pi*i*passo);    % Frequencia de excitaçao (rad/s)
                x = zeros(glf,nt+1);    % Inicializaçao da matriz de deslocamento
                v = zeros(glf,nt+1);    % Inicializaçao da matriz de velocidade
                a = zeros(glf,nt+1);    % Inicializaçao da matriz de aceleraçao
                dx = zeros(glf,1);      % Incremento da matriz de deslocamento
                dv = zeros(glf,1);      % Incremento da matriz de velocidade
                da = zeros(glf,1);      % Incremento da matriz de aceleração
                v(:,1) = zeros(glf,1);  % Velocidade inicial: {v}o = 0
                x(:,1) = zeros(glf,1);  % Deslocamento inicial: {x}o = 0
                Kef = zeros(glf,glf);   % Matriz de rigidez efetiva
                dFef = zeros(glf,nt);   % Vetor de força efetiva
                P = zeros(glf,nt);      % Inicializaçao da vetor força
                Po = 1;                 % Amplitude do carregamento 
            %
            % CARREGAMENTO APLICADO
                % Força aplicada na estrutura
                    if tc == 1            
                        for j = 1:gle               
                          for it = 0:1:nt
                              P(j,it+1) = Po*sin(we*(it)*dt);
                          end     
                        end
                %
                % Força aplicada na base da estrutura 
                    elseif tc == 2
                        U = ones(glf,1);
                        for j = 1:Ntlcd
                            U(gle+j,1)= 0;
                        end
                        for it=0:1:nt    
                            P(:,it+1) = -Mf*U*(sin(we*(it)*dt));
                        end
                    end  
                %
                % Metodo de Newmark
                % Entradas para iniciar o algoritmo
                    a(:,1) = Mf\(P(:,1) - Cf*v(:,1) - Kf*x(:,1));                
                    for it = 1:1:nt
                        if it == 1
                            for j = 1:Ntlcd
                                Cf(gle+j,gle+j) = 0;
                            end
                        else
                            for j = 1:Ntlcd
                                Re = (D/visc)*abs(v(gle+j,it));
                                c = ((rug/(3.7*D))^1.11) - (5.16/Re)*log10((rug/(3.7*D))+(5.09/(Re^0.87)));
                                if c < 0 
                                    f = 0;
                                else
                                    f = (1/(-2*log10(c)))^2;
                                end
                                Cf(gle+j,gle+j) = (ck*f + cdk)*abs(v(gle+j,it)); 
                            end    
                        end  
                        Kef =  a1*Mf + a2*Cf + Kf;
                        dFef = (P(:,it+1)-P(:,it)) + (a3*Mf+ a4*Cf)*v(:,it) + (a5*Mf+ a6*Cf)*a(:,it);
                        dx = Kef\dFef;                    
                        dv = a2*dx - a4*v(:,it) - a6*a(:,it);
                        da = a1*dx - a3*v(:,it) - a5*a(:,it);
                        x(:,it+1) = x(:,it) + dx;
                        v(:,it+1) = v(:,it) + dv;
                        a(:,it+1) = a(:,it) + da; 
                    end 
                    % Deslocamento máximo
%                     for l = 1:gle
                        if abs(max(x(gle,:))) > abs(min(x(gle,:)))
                            xm(i+1,gle) = max(x(gle,:));
                        else
                            xm(i+1,gle) = abs(min(x(gle,:)));
                        end
%                     end
                    % Velocidade máxima
%                     for l = 1:glf
%                         if abs(max(v(l,:))) > abs(min(v(l,:)))
%                             vm(i+1,l) = max(v(l,:));
%                         else
%                             vm(i+1,l) = abs(min(v(l,:)));
%                         end                
%                     end
                    % Aceleração máxima
%                     for l = 1:glf
%                         if abs(max(a(l,:))) > abs(min(a(l,:)))
%                             am(i+1,l) = max(a(l,:));
%                         else
%                             am(i+1,l) = abs(min(a(l,:)));
%                         end                  
%                     end
        end
    end

    %      
    % RESPOSTAS: GRAFICOS 
        if tp == 1 % Espectro de resposta em termos de deslocamento     
            plot(faf(:,1),xm(:,gle))
            xlabel('f (Hz)')
            ylabel('|x| (m)')
            grid
        elseif tp == 2  % Espectro de resposta em termos de velocidade 
            plot(faf(:,1),vm(:,:))
            xlabel('FREQUÊNCIA CARREGAMENTO (Hz)')
            ylabel('VELOCIDADE ABSOLUTA (m/s)')
            grid
        elseif tp == 3  % Espectro de respostas em termos de aceleração  
            plot(faf(:,1),am(:,:))
            xlabel('FREQUÊNCIA CARREGAMENTO (Hz)')
            ylabel('ACELERAÇÃO ABSOLUTA (m/s²)')
            grid    
        end  
     box('on');
     ax = gca;
     ay = gca;
     ax.XAxis.MinorTick = 'on';
     ay.YAxis.MinorTick = 'on';
     grid on;
     ax.XMinorGrid = 'on';
     ay.YMinorGrid = 'on';
%        
elseif ta == 2      
%% RESPOSTA NO DOMINIO DO TEMPO
%
% PASSO NO TEMPO(para tc=1 e tc =2)
    dt = 0.01;  % Passo
    ti = 0;            % Inicio
    tf = 40;           % Fim
%
% ACELERAÇÃO DOS SISMOS
    if tc == 3
        if sismo == 1
            'El Centro' %#ok<NOPTS>
            arquivo = fopen('EL CENTRO.txt','rt');
            [sismo, contador] = fscanf(arquivo, '%f');
            fclose(arquivo);  
            dt = 0.02;
            tf = 35;
        elseif sismo == 2
            'Northridge'    %#ok<NOPTS>
            arquivo2 = fopen('NORTHRIDGEML.txt','rt');
            [sismo,contador2] = fscanf(arquivo2, '%f');
            fclose(arquivo2);  
            dt = 0.005; 
            tf = 40;
        elseif sismo == 3   
            %'Kobe'  %#ok<NOPTS>
            arquivo3 = fopen('KOBETAKATORI90.txt','rt');
            [sismo, contador3] = fscanf(arquivo3, '%f');
            fclose(arquivo3);  
            dt = 0.01;
            tf = 40;
        end
        y = zeros(length(sismo)+1,1);
        t = zeros(length(sismo)+1,1);
        for it=1:1:length(sismo)    
            y(it+1,1) = 9.806*sismo(it,1);
            t(it+1,1) = dt*it;
        end
    end
%
% INTERVALO DA ANÁLISE NO DOMINIO DO TEMPO 
    nt = (tf-ti)/dt;    % Quantidade de passos
    we = 6.0;           % Frequencia de excitaçao (rad/s)    
%
% PARAMETROS CONSTANTES PARA A INTEG NUM PELO MET DE NEWMARK
    gamma = 0.5;
    beta = 0.25;
    a1 = 1/(beta*(dt^2));
    a2 = gamma/(beta*dt);
    a3 = 1/(beta*dt);
    a4 = gamma/beta;
    a5 = 1/(2*beta);
    a6 = (dt)*((gamma/(2*beta))-1);    
%    
% MODELO MATEMÁTICO
    if tm == 1    % Modelo de Shear Building engastado
        % 
        % ALGORITMO: METODO DE NEWMARK (ACELERAÇÃO MEDIA - CHOPRA P. 614) 
        a = zeros(gle,nt+1);    % Inicializaçao da matriz de aceleraçao
        da = zeros(gle,1);      % Incremento da matriz de aceleração
        v = zeros(gle,nt+1);    % Inicializaçao da matriz de velocidade
        dv = zeros(gle,1);      % Incremento da matriz de velocidade
        x = zeros(gle,nt+1);    % Inicializaçao da matriz de deslocamento
        dx = zeros(gle,1);      % Incremento da matriz de deslocamento
        v(:,1) = zeros(gle,1);  % Velocidade inicial: {v}o = 0
        x(:,1) = zeros(gle,1);  % Deslocamento inicial: {x}o = 0
        P = zeros(gle,nt+1);    % Inicializaçao da vetor força 
        U = ones(gle,1);        % Vetor unitário
        dFef = zeros(gle,1);    % Incremento do vetor força
        Po = 1;                 % Amplitude do carregamento 
        %
        % CARREGAMENTO APLICADO NA ESTRUTURA {P}         
            if tc == 1
              for i = 1:gle
                for j = 0:1:nt
                    P(i,j+1) = Po*sin(we*(j)*dt);
                end
              end
            %  
            elseif tc == 2
              U=ones(gle,1);
              for it=0:1:nt    
                  P(:,it+1) = Me*U*(Po*sin(we*(it)*dt));
              end
            %
            elseif tc == 3   
                for it=1:1:length(sismo)    
                    P(:,it) = -Me*U*(9.806)*sismo(it,1);
                end
            end        
        %
        % Metodo de Newmark
        % Entradas para iniciar o algoritmo
            a(:,1) = Me\(P(:,1) - Ce*v(:,1) - Ke*x(:,1));
            Kef =  a1*Me + a2*Ce+ Ke;
            for i = 1:1:nt
                if i*dt <= (tf + dt)
                    dFef = (P(:,i+1)-P(:,i)) + (a3*Me+ a4*Ce)*v(:,i) + (a5*Me+ a6*Ce)*a(:,i) ;
                else
                    P(:,i+1) = zeros(gle,1);
                    dFef = (P(:,i+1)-P(:,i)) + (a3*Me+ a4*Ce)*v(:,i) + (a5*Me+ a6*Ce)*a(:,i) ;
                end
                    dx = Kef\dFef;                    
                    dv = a2*dx - a4*v(:,i) - a6*a(:,i);
                    da = a1*dx - a3*v(:,i) - a5*a(:,i);
                    x(:,i+1) = x(:,i) + dx;
                    v(:,i+1) = v(:,i) + dv;
                    a(:,i+1) = a(:,i) + da; 
            end
            
    % Modelo de Shear Building + ISE
    elseif tm == 2
            x = zeros(glse,nt+1);       % Inicializaçao da matriz de deslocamento
            v = zeros(glse,nt+1);       % Inicializaçao da matriz de velocidade
            a = zeros(glse,nt+1);       % Inicializaçao da matriz de aceleraçao
            dx = zeros(glse,1);         % Incremento da matriz de deslocamento
            dv = zeros(glse,1);         % Incremento da matriz de velocidade
            da = zeros(glse,1);         % Incremento da matriz de aceleração
            u = zeros(glse-2,nt+1);     % Inicializaçao da matriz de deslocamento horizontal total
            at = zeros(glse-2,nt+1);    % Inicializaçao da matriz de aceleraçao total
            vt = zeros(glse-2,nt+1);    % Inicializaçao da matriz de velocidade total
            v(:,1) = zeros(glse,1);     % Velocidade inicial: {v}o = 0
            x(:,1) = zeros(glse,1);     % Deslocamento inicial: {x}o = 0
            dFef = zeros(glse,nt+1);    % Vetor de força efetiva
            P = zeros(glse,nt+1);       % Inicializaçao da vetor força
            Po = 1;                   % Amplitude do carregamento   
            Meff = zeros(glse,1);
            U = zeros(glse,1);
            U(2,1) = 1;
            Ft = zeros(1,nt);
            Mt = zeros(1,nt);          
            %
            % CARREGAMENTO APLICADO NA ESTRUTURA {P} 
            if tc == 1
                for i = 1:gle               
                  for it = 0:1:nt
                      P(i+2,it+1) = Po*sin(we*(it)*dt);
                      P(1,it+1) = Mt(1,it+1) + P(i+2,it+1)*i*leng;
                      P(2,it+1) = Ft(1,it+1) + P(i+2,it+1);
                  end
                  Mt = P(1,:);  
                  Ft = P(2,:);     
                end
            %
            elseif tc == 2                  
                Meff = Mse*U;
                for it=0:1:nt    
                    P(:,it+1) = Meff*(sin(we*(it)*dt));
                end
            %   
            elseif tc == 3                               
                Meff = Mse*U;    
                for it=1:1:length(sismo)    
                    P(:,it) = -Meff*(9.806)*sismo(it,1);
                end       
            end
            %
            % Metodo de Newmark
            % Entradas para iniciar o algoritmo           
                a(:,1) = Mse\(P(:,1) - Cse*v(:,1) - Kse*x(:,1));
                Kef =  a1*Mse + a2*Cse+ Kse;
                for i = 1:1:nt
                    if i*dt <= (tf + dt)
                        dFef = (P(:,i+1)-P(:,i)) + (a3*Mse+ a4*Cse)*v(:,i) + (a5*Mse+ a6*Cse)*a(:,i) ;
                    else
                        P(:,i+1) = zeros(glse,1);
                        dFef = (P(:,i+1)-P(:,i)) + (a3*Mse+ a4*Cse)*v(:,i) + (a5*Mse+ a6*Cse)*a(:,i) ;
                    end
                        dx = Kef\dFef;                    
                        dv = a2*dx - a4*v(:,i) - a6*a(:,i);
                        da = a1*dx - a3*v(:,i) - a5*a(:,i);
                        x(:,i+1) = x(:,i) + dx;
                        v(:,i+1) = v(:,i) + dv;
                        a(:,i+1) = a(:,i) + da;             
                %
                % Resposta do deslocamento horizontal total do sistema acoplado
                   %for j = 1:gle
                   %    u(j,:)= x(j+2,:) + x(2,:) + j*leng*x(1,:);
                   %end
                % Resposta absoluta em termos de velocidade do sistema acoplado
%                    for j = 1:gle
%                        vt(j,:)= v(j+2,:) + v(2,:) + j*leng*v(1,:);
%                    end
                % Resposta absoluta em termos de velocidade do sistema acoplado
                   %for j = 1:gle
                   %    at(j,:)= a(j+2,:) + a(2,:) + j*leng*a(1,:);
                   %end
                end
    % Modelo de Shear Building + TLCD
    elseif tm == 3               
            x = zeros(glf,nt+1);     % Inicializaçao da matriz de deslocamento
            v = zeros(glf,nt+1);     % Inicializaçao da matriz de velocidade
            a = zeros(glf,nt+1);     % Inicializaçao da matriz de aceleraçao
            dx = zeros(glf,1);       % Incremento da matriz de deslocamento
            dv = zeros(glf,1);       % Incremento da matriz de velocidade
            da = zeros(glf,1);       % Incremento da matriz de aceleração
            v(:,1) = zeros(glf,1);   % Velocidade inicial: {v}o = 0
            x(:,1) = zeros(glf,1);   % Deslocamento inicial: {x}o = 0
            Kef = zeros(glf,glf);    % Matriz de rigidez efetiva
            dFef = zeros(glf,nt);    % Vetor de força efetiva
            P = zeros(glf,nt);       % Inicializaçao da vetor força
            Po = 0.2;                % Amplitude do carregamento  
            U = ones(glf,1);         % Distribuição do carregamento
            for i = 1:Ntlcd
                U(gle+i,1)= 0;
            end
        %    
        % CARREGAMENTO APLICADO NA ESTRUTURA {P} 
            if tc == 1
              for i = 1:gle
                for j = 0:1:nt
                    P(i,j+1) = Po*sin(we*(j)*dt);
                end
              end
            %  
            elseif tc == 2
              for it=0:1:nt    
                  P(:,it+1) = -Mf*U*(Po*sin(we*(it)*dt));
              end
            %
            elseif tc == 3 
                for it=1:1:length(sismo)    
                    P(:,it) = -Mf*U*(9.806)*sismo(it,1);
                end
            end        
            %
            % Metodo de Newmark
            % Entradas para iniciar o algoritmo
                a(:,1) = Mf\(P(:,1) - Cf*v(:,1) - Kf*x(:,1));                
                for i = 1:1:nt-1
                    if i == 1
                        for j = 1:Ntlcd
                            Cf(gle+j,gle+j) = 0;
                        end
                    else
                        for j = 1:Ntlcd
                            Re = (D/visc)*abs(v(gle+j,i));
                            c = ((rug/(3.7*D))^1.11) - (5.16/Re)*log10((rug/(3.7*D))+(5.09/(Re^0.87)));
                            if c < 0 
                                f = 0;
                            else
                                f = (1/(-2*log10(c)))^2;
                            end
                            Cf(gle+j,gle+j) = (ck*f + cdk)*abs(v(gle+j,i)); 
                        end
                    end  
                    Kef =  a1*Mf + a2*Cf + Kf;
                    if i*dt <= (tf + dt)
                        dFef = (P(:,i+1)-P(:,i)) + (a3*Mf+ a4*Cf)*v(:,i) + (a5*Mf+ a6*Cf)*a(:,i) ;
                    else
                        P(:,i+1) = zeros(glf,1);
                        dFef = (P(:,i+1)-P(:,i)) + (a3*Mf+ a4*Cf)*v(:,i) + (a5*Mf+ a6*Cf)*a(:,i) ;
                    end
                    dx = Kef\dFef;
                    dv = a2*dx - a4*v(:,i) - a6*a(:,i);
                    da = a1*dx - a3*v(:,i) - a5*a(:,i);
                    x(:,i+1) = x(:,i) + dx;     
                    v(:,i+1) = v(:,i) + dv;
                    a(:,i+1) = a(:,i) + da; 
                end
                
    % Modelo de ISE + Shear Building + TLCD
    elseif tm == 4               
        a = zeros(glia,nt+1);       % Inicializaçao da matriz de aceleraçao
        da = zeros(glia,1);         % Incremento da matriz de aceleração
        v = zeros(glia,nt+1);       % Inicializaçao da matriz de velocidade
        dv = zeros(glia,1);         % Incremento da matriz de velocidade
        x = zeros(glia,nt+1);       % Inicializaçao da matriz de deslocamento
        dx = zeros(glia,1);         % Incremento da matriz de deslocamento
        v(:,1) = zeros(glia,1);     % Velocidade inicial: {v}o = 0
        x(:,1) = zeros(glia,1);     % Deslocamento inicial: {x}o = 0
        Kef = zeros(glia,glia);     % Matriz de rigidez efetiva
        dFef = zeros(glia,nt);      % Vetor de força efetiva
        P = zeros(glia,nt);         % Inicializaçao da vetor força
        Po = 0.2;                   % Amplitude do carregamento  
        Meff = zeros(glia,1);
        U = zeros(glia,1);
        U(2) = 1;
        Ft = zeros(1,nt);
        Mt = zeros(1,nt);          
        %
        % CARREGAMENTO APLICADO NA ESTRUTURA {P} 
            if tc == 1
                for i = 1:gle               
                  for it = 0:1:nt
                      P(i+2,it+1) = Po*sin(we*(it)*dt);
                      P(1,it+1) = Mt(1,it+1) + P(i+2,it+1)*i*leng;
                      P(2,it+1) = Ft(1,it+1) + P(i+2,it+1);
                  end
                  Mt = P(1,:);  
                  Ft = P(2,:);     
                end
            %
            elseif tc == 2   
                Meff = Mia*U;
                for it=0:1:nt    
                    P(:,it+1) = -Po*Meff*(sin(we*(it)*dt));
                end
            %   
            elseif tc == 3 
                Meff = Mia*U;   
                for it=1:1:length(sismo)    
                    P(:,it) = -Meff*(9.806)*sismo(it,1);
                end       
            end
            %
            % Metodo de Newmark
            % Entradas para iniciar o algoritmo
                a(:,1) = Mia\(P(:,1) - Cia*v(:,1) - Kia*x(:,1));                
                for i = 1:1:nt
                    if i == 1
                        for j = 1:Ntlcd
                            Cia(glse+j,glse+j) = 0;
                        end
                    else
                        for j = 1:Ntlcd
                            Re = (D/visc)*abs(v(glse+j,i));
                            c = ((rug/(3.7*D))^1.11) - (5.16/Re)*log10((rug/(3.7*D))+(5.09/(Re^0.87)));
                            if c < 0 
                                f = 0;
                            else
                                f = (1/(-2*log10(c)))^2;
                            end
                            Cia(glse+j,glse+j) = (ck*f + cdk)*abs(v(glse+j,i));
                        end
                    end  
                    Kef =  a1*Mia + a2*Cia + Kia;
                    dFef = (P(:,i+1)-P(:,i)) + (a3*Mia+ a4*Cia)*v(:,i) + (a5*Mia+ a6*Cia)*a(:,i);
                    dx = Kef\dFef;
                    dv = a2*dx - a4*v(:,i) - a6*a(:,i);
                    da = a1*dx - a3*v(:,i)- a5*a(:,i);
                    x(:,i+1) = x(:,i) + dx;     
                    v(:,i+1) = v(:,i) + dv;
                    a(:,i+1) = a(:,i) + da; 
                end            
    end
 %
 %% GRAFICO: RESPOSTA DINAMICA 
%      if tc == 3
%         subplot(211)
%         plot(t(1:nt,1),y(1:nt,1))
%         xlabel('t(s)')
%         ylabel('Xg (m/s²)')      
%         if sismo == 1
%             title('Acelerograma - El Centro')
%         elseif sismo == 2
%             title('Acelerograma - Northridge')
%         elseif sismo == 3
%             title('Acelerograma - Kobe')
%         elseif sismo == 4
%             title('Acelerograma - Tabas')
%         end
%      elseif tc == 1
%         time = 0:dt:(nt*dt);
%         subplot(211)
%         plot(time,P) 
%         xlabel('t(s)')
%         ylabel('f(t)')
%         title('FORÇA DINÂMICA APLICADA NA ESTRUTURA')
%       elseif tc == 2
%         time = 0:dt:(nt*dt);
%         subplot(211)
%         plot(time,P) 
%         xlabel('t(s)')
%         ylabel('f(t)')
%         title('ACELERAÇÃO APLICADA NA BASE ESTRUTURA')
%      end
%      box('on');
%      ax = gca;
%      ay = gca;
%      ax.XAxis.MinorTick = 'on';
%      ay.YAxis.MinorTick = 'on';
%      grid on;
%      ax.XMinorGrid = 'on';
%      ay.YMinorGrid = 'on';
%    
% Resposta horizontal total do sistema acoplado
        if tm == 1
            n = gle;
        elseif tm == 2
            n = glse;
        elseif tm == 3
            n = glf-Ntlcd;
        elseif tm == 4
            n = glia-Ntlcd;
        end
        time = 0:dt:(nt*dt);
        [pks,loc] = findpeaks(x(n,:),'MinPeakDistance',50);
        if tp == 1 
        % Plotagem da resposta em termos de deslocamento
%             subplot(212)
            plot(time,x(n,:),'black',loc*0.01,pks,'o')                       
            %title('Vibração forçada')
            xlabel('t (s)')
            ylabel('x (m)')
            grid 

        elseif tp == 2
        % Plotagem da resposta em termos de velocidade
            subplot(212)
            plot(time,v(n,:))
            xlabel('t (s)')
            ylabel('VELOCIDADE (m/s)')
            grid

        elseif tp == 3
        % Plotagem da resposta em termos de aceleração;
            subplot(212)
            plot(time,a(n,:))
            xlabel('t (s)')
            ylabel('ACELERAÇÃO (m/s²)')
            grid
        end
        ax = gca;
        ay = gca;
        ax.XAxis.MinorTick = 'on';
        ay.YAxis.MinorTick = 'on';
        grid on;
        ax.XMinorGrid = 'on';
        ay.YMinorGrid = 'on';
end

%% RESULTADOS
%
% FREQUENCIA
    %    fse(1)
    %    1/fse(1)
    %    fn(1)/fse(1)-1

% PARAMETROS ATENUADOR DE COLUNA LÍQUIDA
    if tm == 3
        %wn(1:3)                                                     %#ok<NOPTS>       
        %ftlcd                                                      %#ok<NOPTS>
        %wtlcd                                                      %#ok<NOPTS> 
        wptlcd                                                      %#ok<NOPTS> 
        Razao_de_Massa_Atenuador_Estrutura = (Ntlcd*Mtlcd/Met)*100 %#ok<NOPTS>
        %Deslocamento_Coluna_Liquida = max(x(n+Ntlcd,:))            %#ok<NOPTS>
        %Deslocamento_Maximo_Permitido = (L-B)/2                    %#ok<NOPTS>                                %#ok<NOPTS>
    end    
    if tm == 4
        %fn(1)                                                       %#ok<NOPTS>
        %ftlcd                                                       %#ok<NOPTS>
        %wtlcd                                                       %#ok<NOPTS> 
        wptlcd                                                      %#ok<NOPTS> 
        Razao_de_Massa_Atenuador_Estrutura = (Ntlcd*Mtlcd/Met)*100  %#ok<NOPTS>
        %Deslocamento_Coluna_Liquida = max(x(n+Ntlcd,:))             %#ok<NOPTS>
        %Deslocamento_Maximo_Permitido = Z                           %#ok<NOPTS>                                %#ok<NOPTS>
    end    
% DESLOCAMENTOS   
%         if max(x(n,:)) > abs(min(x(n,:)))
%             Xmax = max(x(n,:))                                  %#ok<NOPTS>
%         else
%             Xmax = abs(min(x(n,:)))                             %#ok<NOPTS>
%         end
        Xrms = rms(x(n,:))                                      %#ok<NOPTS>
%         if tm == 2
%             if max(u(n-2,:)) > abs(min(u(n-2,:)))
%                 Umax = max(u(n-2,:))                            %#ok<NOPTS>
%             else
%                     Umax = abs(min(u(n-2,:)))                   %#ok<NOPTS>
%             end
%             Urms = rms(u(n-2,:))                                %#ok<NOPTS>
%         end
    
% ACELERAÇÃO    
%         if max(a(n,:)) > abs(min(a(n,:)))
%             Amax = max(a(n,:))
%         else
%             Amax = abs(min(a(n,:)))
%         end
%         Arms = rms(a(n,:))    
%         if tm == 2
%             if max(at(n-2,:)) > abs(min(at(n-2,:)))
%                 Atmax = max(at(n-2,:))
%             else
%                 Atmax = abs(min(at(n-2,:)))
%             end
%             Atrms = rms(at(n-2,:))
%         end
    

