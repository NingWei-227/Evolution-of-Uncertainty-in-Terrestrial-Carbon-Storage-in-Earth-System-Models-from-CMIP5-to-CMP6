clear
clc
% Figure 4: Attributing the variance of global terrestrial C storage among models into different components
% global terrestrial carbon storage (X)
% net primary production (NPP)
% Ecosystem carbon residence time (tuaE)
% carbon storage potential (Xp)
% gross primary production (GPP)
% carbon use efficiency (CUE)
% baseline ecosystem carbon residence time (tuaE_base)
% environmental scalar, temperature (tuaE_scalar_tem)
% environmental scalar, precipitation (tuaE_scalar_pre)

% CMIP5: 
% constributions of NPP, tuaE, and Xp to X
% NPP_2X = 0.498; tuaE_2X = 0.502; 
% contributions of GPP and CUE to NPP
% GPP_2NPP = 0.685; CUE_2NPP = 0.315 
% contributions of tuaE_base, tuaE_scalar_tem, tuaE_scalar_pre to tuaE
% tuaE_base = 0.98
% tuaE_scalar_tem = 0.01
% tuaE_scalar_pre = 0.01 

% CMIP6: 
% constributions of NPP, tuaE, and Xp to X
% NPP_2X = 0.255; tuaE_2X = 0.745; 
% contributions of GPP and CUE to NPP
% GPP_2NPP = 0.495; CUE_2NPP = 0.505 
% contributions of tuaE_base, tuaE_scalar_tem, tuaE_scalar_pre to tuaE
% tuaE_base = 0.665
% tuaE_scalar_tem = 0.328
% tuaE_scalar_pre = 0.007 

%% 
% define angle and radius of the circular diagrams
theta1=(180-40)/180*pi;
theta2=(180+65)/180*pi;
theta3 = -65/180*pi;
theta4 = 40/180*pi;
the1=theta1:pi/180:theta2;
the2=theta3:pi/180:theta4;

r=5;

x1=r*cos(the1);
y1=r*sin(the1);
x2=r*cos(the2);
y2=r*sin(the2);

% circular diagrams
figure
set(gcf,'position',[100 100 680,440.8],'color',[1 1 1])
panel = tight_subplot(1,2,[0.08 0.16],[0.08 0.34],[0.03 0.03])
cmip5_panel = panel(1)
axes(cmip5_panel)
hold on 
plot(x1,y1,'color',[0.5 0.5 0.5],'LineWidth',2);
plot(x2,y2,'color',[0.5 0.5 0.5],'LineWidth',2)

k1 = tan(19/180*pi);
k2 = -tan(19/180*pi);
line1_x = 0:-0.01:-5.*cos(19/180*pi);
line2_x = 0:0.01:5*cos(19/180*pi);
line1_y = k1*line1_x;
line2_y = k2*line2_x;

plot(line1_x,line1_y,'color',[0.5 0.5 0.5],'LineWidth',2);
plot(line2_x,line2_y,'color',[0.5 0.5 0.5],'LineWidth',2);
%plot([0 0],[0 -2.3],'color',[0.5 0.5 0.5],'LineWidth',2);

% relative contributions were calculated in the R scripts
text(0,0,'cLand','HorizontalAlignment','center',...
        'FontName','Arial','FontSize',14,...
        'BackgroundColor',[1 1 1],'Margin',6)
       
    
text(3.*cos(19/180*pi),-3.*sin(19/180*pi),'\tau_E','HorizontalAlignment','center',...
        'FontName','Arial','FontSize',12,...
        'BackgroundColor',[1 1 1],'Margin',1,'Rotation',78)   
    text(2.5087,0.1944,'0.502','HorizontalAlignment','center',...
        'FontName','Arial','FontSize',11,...
        'Margin',1,'Rotation',68,'color',[0.89,0.02 0.02])       
            
text(-3.*cos(19/180*pi),-3.*sin(19/180*pi),'NPP','HorizontalAlignment','center',...
        'FontName','Arial','FontSize',12,...
        'BackgroundColor',[1 1 1],'Margin',2,'Rotation',-74)  
    text(-2.55,0.4072,'0.498','HorizontalAlignment','center',...
        'FontName','Arial','FontSize',12,...
        'Margin',1,'Rotation',-75,'color',[0.89,0.02 0.02]) 
    
    
text(-5.*cos(45/180*pi),5.*sin(45/180*pi),'GPP','HorizontalAlignment','center',...
        'FontName','Arial','FontSize',12,...
       'Margin',2,'Rotation',-41) 
   text(-4.67,2.9075,'0.68','HorizontalAlignment','center',...
        'FontName','Arial','FontSize',12,...
       'Margin',2,'Rotation',49,'color',[1,0,0],'FontWeight','bold')   
text(-5.*cos(69/180*pi),-5.*sin(69/180*pi),'CUE','HorizontalAlignment','center',...
        'FontName','Arial','FontSize',12,...
        'BackgroundColor',[1 1 1],'Margin',2,'Rotation',-10) 
    text(-3.07,-4.866,'0.32','HorizontalAlignment','center',...
        'FontName','Arial','FontSize',11,...
       'Margin',2,'Rotation',-30,'color',[0.89,0.02 0.02]) 
       
      
text(5.*cos(19/180*pi),-5.*sin(19/180*pi),'\tau_E','HorizontalAlignment','center',...
        'FontName','Arial','FontSize',12,...
        'BackgroundColor',[1 1 1],'Margin',1,'Rotation',78)
    text(5.41,-0.545,'0.98','HorizontalAlignment','center',...
        'FontName','Arial','FontSize',12,...
        'Margin',1,'Rotation',78,'color',[1,0,0],'FontWeight','bold')
    
 text(5.*cos(45/180*pi),5.*sin(45/180*pi),'\xi_w','HorizontalAlignment','center',...
        'FontName','Arial','FontSize',13,...
       'Margin',2,'Rotation',41)  
   text(4.67,2.9075,'0.01','HorizontalAlignment','center',...
        'FontName','Arial','FontSize',11,...
       'Margin',2,'Rotation',-49,'color',[0.89,0.02 0.02]) 
 text(5.*cos(69/180*pi),-5.*sin(69/180*pi),'\xi_T','HorizontalAlignment','center',...
        'FontName','Arial','FontSize',13,...
        'BackgroundColor',[1 1 1],'Margin',2,'Rotation',10) 
    text(2.7033,-5.202,'0.01','HorizontalAlignment','center',...
        'FontName','Arial','FontSize',11,...
       'Margin',2,'Rotation',30,'color',[0.89,0.02 0.02])   
axis([-5 5 -5 5])
set(gca,'YAxisLocation','origin');% set y axis center around the origin
set(gca,'XAxisLocation','origin');% set x axis center around the origin
set(gca,'xcolor','none','ycolor','none')   


    
       
% CMIP6
cmip6_panel = panel(2)
axes(cmip6_panel)
hold on
plot(x1,y1,'color',[0.5 0.5 0.5],'LineWidth',2);
plot(x2,y2,'color',[0.5 0.5 0.5],'LineWidth',2)

k1 = tan(19/180*pi);
k2 = -tan(19/180*pi);
line1_x = 0:-0.01:-5.*cos(19/180*pi);
line2_x = 0:0.01:5*cos(19/180*pi);
line1_y = k1*line1_x;
line2_y = k2*line2_x;

plot(line1_x,line1_y,'color',[0.5 0.5 0.5],'LineWidth',2);
plot(line2_x,line2_y,'color',[0.5 0.5 0.5],'LineWidth',2);
%plot([0 0],[0 -2.3],'color',[0.5 0.5 0.5],'LineWidth',2);

text(0,0,'cLand','HorizontalAlignment','center',...
        'FontName','Arial','FontSize',14,...
        'BackgroundColor',[1 1 1],'Margin',6)
    
    
text(2.9.*cos(19/180*pi),-2.9.*sin(19/180*pi),'\tau_E','HorizontalAlignment','center',...
        'FontName','Arial','FontSize',12,...
        'BackgroundColor',[1 1 1],'Margin',1,'Rotation',74) 
    text(2.5087,0.1944,'0.74','HorizontalAlignment','center',...
        'FontName','Arial','FontSize',12,...
        'Margin',1,'Rotation',68,'color',[1,0,0],'FontWeight','bold') 
        
text(-2.9.*cos(19/180*pi),-2.9.*sin(19/180*pi),'NPP','HorizontalAlignment','center',...
        'FontName','Arial','FontSize',12,...
        'BackgroundColor',[1 1 1],'Margin',2,'Rotation',-74) 
    text(-2.55,0.4072,'0.26','HorizontalAlignment','center',...
        'FontName','Arial','FontSize',11,...
        'Margin',1,'Rotation',-75,'color',[0.89,0.02 0.02]) 
    
    
text(-5.*cos(45/180*pi),5.*sin(45/180*pi),'GPP','HorizontalAlignment','center',...
        'FontName','Arial','FontSize',12,...
       'Margin',2,'Rotation',-41)
   text(-4.67,2.9075,'0.50','HorizontalAlignment','center',...
        'FontName','Arial','FontSize',11,...
       'Margin',2,'Rotation',49,'color',[0.89,0.02 0.02])     
text(-5.*cos(69/180*pi),-5.*sin(69/180*pi),'CUE','HorizontalAlignment','center',...
        'FontName','Arial','FontSize',12,...
        'BackgroundColor',[1 1 1],'Margin',2,'Rotation',-10) 
    text(-3.07,-4.866,'0.50','HorizontalAlignment','center',...
        'FontName','Arial','FontSize',11,...
       'Margin',2,'Rotation',-30,'color',[0.89,0.02 0.02])
    
           
text(5.*cos(19/180*pi),-5.*sin(19/180*pi),'\tau_E','HorizontalAlignment','center',...
        'FontName','Arial','FontSize',12,...
        'BackgroundColor',[1 1 1],'Margin',1,'Rotation',74)
    text(5.41,-0.545,'0.67','HorizontalAlignment','center',...
        'FontName','Arial','FontSize',12,...
        'Margin',1,'Rotation',78,'color',[1,0,0],'FontWeight','bold')
    
 text(5.*cos(45/180*pi),5.*sin(45/180*pi),'\xi_w','HorizontalAlignment','center',...
        'FontName','Arial','FontSize',13,...
       'Margin',2,'Rotation',41) 
   text(4.67,2.9075,'0.01','HorizontalAlignment','center',...
        'FontName','Arial','FontSize',11,...
       'Margin',2,'Rotation',-49,'color',[0.89,0.02 0.02])
 text(5.*cos(69/180*pi),-5.*sin(69/180*pi),'\xi_T','HorizontalAlignment','center',...
        'FontName','Arial','FontSize',13,...
        'BackgroundColor',[1 1 1],'Margin',2,'Rotation',10)
    text(2.7033,-5.202,'0.32','HorizontalAlignment','center',...
        'FontName','Arial','FontSize',11,...
       'Margin',2,'Rotation',30,'color',[0.89,0.02 0.02])   
axis([-5 5 -5 5])
set(gca,'YAxisLocation','origin');% set y axis center around ordinate origin
set(gca,'XAxisLocation','origin');% set x axis center around ordinate origin
set(gca,'xcolor','none','ycolor','none')



% pie plots 
% CMIP5
pies_panel = tight_subplot(1,2,[0.08 0.16],[0.53 0.04],[0.03 0.03])
NPP_2X = 0.498;
tuaE_2X = 0.502;

GPP = 0.498*0.68
CUE = 0.498*0.32

tuaE_base = 0.502*0.98
tuaE_scalar_tem = 0.502*0.01
tuaE_scalar_pre = 0.502*0.01 

axes(pies_panel(1))
CMIP5_bg = [GPP CUE tuaE_base tuaE_scalar_tem tuaE_scalar_pre]
labels = {'GPP' 'CUE' '\tau_E' '\xi_T' '\xi_w'}
h5_ed5 = pie(CMIP5_bg,labels)
h5_ed5(1).EdgeColor = [1 1 1];
h5_ed5(3).EdgeColor = [1 1 1];
h5_ed5(5).EdgeColor = [1 1 1];
h5_ed5(7).EdgeColor = [1 1 1];
h5_ed5(9).EdgeColor = [1 1 1];

colormap([119 172 48;
          204 204 0;
          0 191 191;
          255 102 102;
          0 0 255]./255)      
text(-1.313,1.02,'(a) CMIP5','HorizontalAlignment','center',...
        'FontName','Arial','FontSize',13)        
    
% pie plots 
% CMIP6    
axes(pies_panel(2))
NPP_2X = 0.26
tuaE_2X = 0.74
GPP = 0.26*0.5
CUE = 0.26*0.5

tuaE_base = 0.74*0.67
tuaE_scalar_tem = 0.74*0.32
tuaE_scalar_pre = 0.74*0.01

CMIP6_bg = [GPP CUE tuaE_base tuaE_scalar_tem tuaE_scalar_pre]
labels = {'GPP' 'CUE' '\tau_E' '\xi_T' '\xi_w'}
h6_ed5 = pie(CMIP6_bg,labels)

colormap([119 172 48;
          204 204 0;
          0 191 191;
          255 102 102;
          0 0 255]./255)      
h6_ed5(1).EdgeColor = [1 1 1];
h6_ed5(3).EdgeColor = [1 1 1];
h6_ed5(5).EdgeColor = [1 1 1];
h6_ed5(7).EdgeColor = [1 1 1];
h6_ed5(9).EdgeColor = [1 1 1];    
text(-1.313,1.02,'(b) CMIP6','HorizontalAlignment','center',...
        'FontName','Arial','FontSize',13)       
    
