%----------------------------------PROYECTO FINAL
   n = 48;
  
%PROPIEDADES DEL MATERIAL (COBRE)
alpha_expansion = 1.67E-5; % Thermal Expansion Coefficient [1/C]
rho = 8940; % Density [kg/m^3]
k = 401; % Thermal Conductivity [W/m C]
C =  450; % Specific Heat [kg C]

%OBTENCION PROPIEDADES DERIVDAS
    %Thermal Difussivity
alpha = k/(rho * C); %[mm^2/s]

% COEFICIENTE DE CONVECCION DEL FLUIDO
    %AIRE (CONVECCION NATURAL)
        %TEMPERATURA ASUMIDA DE 40ºC
H_out = 100 ; % [W/m^2K]
T_aire_out = 40 ; % [C]

    %AIRE (CONVECCION NATURAL)
        %TEMPERATURA ASUMIDA DE 50ºC
H_in = 200 ; % [W/m^2K]
T_aire_in = 50; % [C] 

%GEOMETRÍA
    %LONGITUD DE LA MATRIZ
        X= 6; %[m]
        Y= 6; %[m]
    
    %LONGITUD RECTANGULO GRANDE
    
    %LONGITUD RECTANGULO
        x = 2; %[m] cuadrado pequeño 
        y = 1; %[m] cuadrado pequeño 
        
    %DISTANCIAS ENTRE NODOS
        dx = X/(n-1);
        dy = Y/(n-1);

%NUMERO DE BIOT OUT
Bi_out = (H_out*dx)/k;

%NUMERO DE BIOT IN
Bi_in = (H_in*dx)/k;

%OBTENCION DELTA t 

    %FORMULA 1 (PARA FOURIER <= 1/4)
        dt1 = (1/4)* (dx^2)/alpha;
        
 %OUT 
    % FORMULA 2(PARA FOURIER (3+BIOT) <= 3/4)
        dt2_out = (3/4)*(dx^2)/(alpha*(3+((H_out*dx)/k)));
        
    % FORMULA 3(PARA FOURIER (2+BIOT) <= 1/2)
        dt3_out = (1/2)*(dx^2)/(alpha*(2+(H_out*dx)/k));

 %IN
    
    % FORMULA 2(PARA FOURIER (3+BIOT) <= 3/4)
        dt2_in = (3/4)*(dx^2)/(alpha*(3+((H_in*dx)/k)));
        
    % FORMULA 3(PARA FOURIER (2+BIOT) <= 1/2)
        dt3_in = (1/2)*(dx^2)/(alpha*(2+(H_in*dx)/k));
        
        V_min = [dt1,dt2_in,dt2_out,dt3_in,dt3_out];
        
        dt = 0.95*min(V_min);
        
%NUMERO DE FOURIER

    %PARA FORMULA 1
        Fo = (alpha*dt)/(dx^2);
        

%----------------------------------CREACION MATRIZ
p = 5000;
T = zeros(n,n,p);

% CONDICIONES INICIALES
% Condicion inicial de 0 en la primera fila

for i = 1:p
    T(1,n/2-8:n/2+8,i) = 0; 
end

%Condicion inicial para los extremos izquierdos, derechos e inferiores de
%40 de temperatura

for i = 1:p
    
    T(9:48,1,i) = 25; %Llenar linea de la izquierda
    T(9:48,48,i) = 25; %Llenar linea dercha
    T(48,2:47,i) = 5; %Llenar linea inferior
    
end

% Condicion inicial de conveccion en las esquinas superiores
    % 1/6*48 = 8 por eso son las 8 primeras filas de la matriz
for i = 1:p
    T(1:8,1:n/2-8,i) = 40;
    T(1:8,n/2+9:n,i) = 40;
end

% Condicion inicial de conveccion en la parte media de la placa
    % 5/6*48 = 40
for i = 1:p
    T(25:32,n/2-7:n/2+8,i) = 50;
end    

%----------------------------------LLENAR MATRIZ EN EL TIEMPO
for i=2:p
%PRIMERA FORMULA (SOLO CONDUCCION - INTERIOR NODE)
    %T(m,n,p)=Fo*(T(m+1,n,p-1)+T(m-1,n,p-1)+T(m,n+1,p-1)+T(m,n-1,p-1))+(1-4*Fo)*T(m,n,p-1)
        %BLOQUE A
             T(2:9,18:31,i) = Fo*(T(3:10,18:31,i-1)+T(1:8,18:31,i-1) + T(2:9,19:32,i-1) + T(2:9,17:30,i-1)) + (1-4*Fo)* T(2:9,18:31,i-1); %es -6 y +6 porque el otro nodo es conveccion
        %BLOQUE B
             T(10:23,2:47,i) = Fo*(T(11:24,2:47,i-1)+T(9:22,2:47,i-1) + T(10:23,3:48,i-1) + T(10:23,1:46,i-1)) + (1-4*Fo)* T(10:23,2:47,i-1);
        %BLQOUE C
             T(24:33,2:15,i) = Fo*(T(25:34,2:15,i-1)+T(23:32,2:15,i-1) + T(24:33,3:16,i-1) + T(24:33,1:14,i-1)) + (1-4*Fo)* T(24:33,2:15,i-1);
        %BLQOUE D
             T(24:33,34:47,i) = Fo*(T(25:34,34:47,i-1)+T(23:32,34:47,i-1) + T(24:33,35:48,i-1) + T(24:33,33:46,i-1)) + (1-4*Fo)* T(24:33,34:47,i-1);
        %BLQOUE E
             T(34:47,2:47,i) = Fo*(T(35:48,2:47,i-1)+T(33:46,2:47,i-1) + T(34:47,3:48,i-1) + T(34:47,1:46,i-1)) + (1-4*Fo)* T(34:47,2:47,i-1);
             
%SEGUNDA FORMULA (CONVECCION EN UNA ESQUINA - NODE AT INTERIOR CORNER WITH CONVECTION)
    %T(m,n,p)=(2/3)*Fo*(2*T(m+1,n,p-1)+2*T(m,n+1,p-1)+T(m,n-1,p-1)+T(m-1,n,p-1)+2*Bi*Tamb)+(1-4*Fo-(4/3)*Bi*Fo)*T(m,n,p-1)
        %PUNTO a
             T(9,17,i) = (2/3)*Fo*(2*T(10,17,i-1)+2*T(9,18,i-1)+T(9,16,i-1)+T(8,17,i-1)+2*Bi_out*T_aire_out)+(1-4*Fo-(4/3)*Bi_out*Fo)*T(9,17,i-1);
        %PUNTO f
             T(33,33,i) = (2/3)*Fo*(2*T(34,33,i-1)+2*T(33,34,i-1)+T(33,32,i-1)+T(32,33,i-1)+2*Bi_in*T_aire_in)+(1-4*Fo-(4/3)*Bi_in*Fo)*T(33,33,i-1);
             
     %T(m,n,p)=(2/3)*Fo*(2*T(m,n-1,p-1)+2*T(m+1,n,p-1)+T(m-1,n,p-1)+T(m,n+1,p-1)+2*Bi*Tamb)+(1-4*Fo-(4/3)*Bi*Fo)*T(m,n,p-1)
         %PUNTO b
             T(9,32,i) = (2/3)*Fo*(2*T(9,31,i-1)+2*T(10,32,i-1)+T(8,32,i-1)+T(9,33,i-1)+2*Bi_out*T_aire_out)+(1-4*Fo-(4/3)*Bi_out*Fo)*T(9,32,i-1);  
     
         %PUNTO e
             T(33,16,i) = (2/3)*Fo*(2*T(33,15,i-1)+2*T(34,16,i-1)+T(32,16,i-1)+T(33,17,i-1)+2*Bi_in*T_aire_in)+(1-4*Fo-(4/3)*Bi_in*Fo)*T(33,16,i-1);
             
      %T(m,n,p)=(2/3)*Fo*(2*T(m-1,n,p-1)+2*T(m,n-1,p-1)+T(m+1,n,p-1)+T(m,n+1,p-1)+2*Bi*Tamb)+(1-4*Fo-(4/3)*Bi*Fo)*T(m,n,p-1)
        %PUNTO c
             T(24,16,i) = (2/3)*Fo*(2*T(23,16,i-1)+2*T(24,15,i-1)+T(25,16,i-1)+T(24,17,i-1)+2*Bi_in*T_aire_in)+(1-4*Fo-(4/3)*Bi_in*Fo)*T(24,16,i-1);
      
      %T(m,n,p)=(2/3)*Fo*(2*T(m-1,n,p-1)+2*T(m,n+1,p-1)+T(m,n-1,p-1)+T(m+1,n,p-1)+2*Bi*Tamb)+(1-4*Fo-(4/3)*Bi*Fo)*T(m,n,p-1)
        %PUNTO d
             T(24,33,i) = (2/3)*Fo*(2*T(23,33,i-1)+2*T(24,34,i-1)+T(24,32,i-1)+T(25,33,i-1)+2*Bi_in*T_aire_in)+(1-4*Fo-(4/3)*Bi_in*Fo)*T(24,33,i-1);
             
        
     
%TERCERA FORMULA (NODE AT PLANE SURFACE WITH CONVECTION)
    %T(m,n,p)=Fo*(2*T(m+1,n,p-1)+T(m,n+1,p-1)+T(m,n-1,p-1)+2*Bi*Tamb)+(1-4*Fo-2*Bi*Fo)*T(m,n,p-1)
        %FILA G
            T(9,2:16,i) = Fo*(2*T(10,2:16,i-1)+T(9,3:17,i-1)+T(9,1:15,i-1)+2*Bi_out*T_aire_out)+(1-4*Fo-2*Bi_out*Fo)*T(9,2:16,i-1);
        %FILA J
            T (9,33:47,i)= Fo*(2*T(10,33:47,i-1)+T(9,34:48,i-1)+T(9,32:46,i-1)+2*Bi_out*T_aire_out)+(1-4*Fo-2*Bi_out*Fo)*T(9,33:47,i-1);
        %FILA N
             T (33,17:32,i) = Fo*(2*T(34,17:32,i-1)+T(33,18:33,i-1)+T(33,16:31,i-1)+2*Bi_in*T_aire_in)+(1-4*Fo-2*Bi_in*Fo)*T(33,17:32,i-1);
           
%     %T(m,n,p)=Fo*(2*T(m,n+1,p-1)+T(m-1,n,p-1)+T(m+1,n,p-1)+2*Bi*Tamb)+(1-4*Fo-2*Bi*Fo)*T(m,n,p-1)        
         %FILA H
             T(2:8,17,i) = Fo*(2*T(2:8,18,i-1)+T(1:7,17,i-1)+T(3:9,17,i-1)+2*Bi_out*T_aire_out)+(1-4*Fo-2*Bi_out*Fo)*T(2:8,17,i-1);
%         %FILA M
            T (25:32,33,i)= Fo*(2*T(25:32,34,i-1)+T(24:31,33,i-1)+T(26:33,33,i-1)+2*Bi_in*T_aire_in)+(1-4*Fo-2*Bi_in*Fo)*T(25:32,33,i-1);
             
%     %T(m,n,p)=Fo*(2*T(m,n-1,p-1)+T(m-1,n,p-1)+T(m+1,n,p-1)+2*Bi*Tamb)+(1-4*Fo-2*Bi*Fo)*T(m,n,p-1)                
        %FILA I
             T (2:8,32,i) = Fo*(2*T(2:8,31,i-1)+T(1:7,32,i-1)+T(3:9,32,i-1)+2*Bi_out*T_aire_out)+(1-4*Fo-2*Bi_out*Fo)*T(2:8,32,i-1);
%       %FILA L
             T (25:32,16,i)= Fo*(2*T(25:32,15,i-1)+T(24:31,16,i-1)+T(26:33,16,i-1)+2*Bi_in*T_aire_in)+(1-4*Fo-2*Bi_in*Fo)*T(25:32,16,i-1);
             
%     %T(m,n,p)=Fo*(2*T(m-1,n,p-1)+T(m,n+1,p-1)+T(m,n-1,p-1)+2*Bi*Tamb)+(1-4*Fo-2*Bi*Fo)*T(m,n,p-1) 
        %FILA K
             T(24,17:32,i) = Fo*(2*T(23,17:32,i-1)+T(24,18:33,i-1)+T(24,16:31,i-1)+2*Bi_in*T_aire_in)+(1-4*Fo-2*Bi_in*Fo)*T(24,17:32,i-1);
         

       
%ERROR PARA CADA BLOQUE DE LA FIGURA
error_bloqueA = norm (T(2:9,18:31,i) - T(2:9,18:31,i-1));
error_bloqueB = norm (T(10:23,2:47,i) -  T(10:23,2:47,i-1));
error_bloqueC = norm (T(24:33,2:15,i) - T(24:33,2:15,i-1));
error_bloqueD = norm (T(24:33,34:47,i) - T(24:33,34:47,i-1));
error_bloqueE = norm (T(34:47,2:47,i) - T(34:47,2:47,i-1));

      %CRITERIO DE CONVERGENCIA
          if (error_bloqueA < 0.0001 && error_bloqueB < 0.0001 && error_bloqueC < 0.0001 && error_bloqueD < 0.0001 && error_bloqueE < 0.0001)
              i
              P_FINAL = i;
              break;
          end  



end

%CODIGO PARA DAR LA VUELTA
T_final = flipud(T);

%SIMULAR
for i=1:P_FINAL
contourf (T_final(:,:, i),12,'edgecolor','none')
colormap(jet) %CHANGE COLOR
axis equal
title ('         Simulación de Transferencia de Calor \newline              por Conducción y Convección   \newline                      en una Placa de Cobre','Fontsize',13);
xlabel ('Nodos', 'FontSize',12);
ylabel ('Nodos', 'FontSize',12);
C = colorbar;   %UBICACION DE LA ESCALA
C.Label.String = 'Temperatura [^{\circ}C]';
caxis([0 25]);
barra = colorbar;
set(barra,'YTick',[0:2:25]); %SACAR EL PASO DE 2 EN 2
drawnow limitrate

end

% %PLOTEAR
% contourf (T_final(:,:, P_FINAL),12,'edgecolor','none')
% colormap(jet(260)) %CHANGE COLOR
% axis equal
% title ('         Simulación de Transferencia de Calor \newline              por Conducción y Convección   \newline                      en una Placa de Cobre','Fontsize',13);
% xlabel ('Nodos', 'FontSize',12);
% ylabel ('Nodos', 'FontSize',12);
% C = colorbar;   %UBICACION DE LA ESCALA
% C.Label.String = 'Temperatura [^{\circ}C]';
% caxis([0 25]);
% barra = colorbar;
% set(barra,'YTick',[0:2:25]); %SACAR EL PASO DE 2 EN 2