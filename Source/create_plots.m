% create plot 
limit_x_start =0;
limit_x_end = 50;
 

f1=figure('units','inch','position',[0,0,8,8]);
subplot(2,2,1); ylabel('y');    xlabel('time');hold on;set(get(gca,'YLabel'),'Rotation',0)
subplot(2,2,2); ylabel('\Pi');  xlabel('time');hold on;set(get(gca,'YLabel'),'Rotation',0)
subplot(2,2,3); ylabel('R');    xlabel('time');hold on;set(get(gca,'YLabel'),'Rotation',0)
subplot(2,2,4); ylabel('\beta');xlabel('time');hold on;set(get(gca,'YLabel'),'Rotation',0)

subplot(2,2,1); plot(0:param.TTT-1,[X_PF(:,1)],'linewidth',2,'color','k','linestyle',':');xlim([limit_x_start,limit_x_end])
subplot(2,2,2); plot(0:param.TTT-1,(ss.Pi+[X_PF(:,2)]).^4-1,'linewidth',2,'color','k','linestyle',':');xlim([limit_x_start,limit_x_end])
subplot(2,2,3); plot(0:param.TTT-1,(ss.R+[X_PF(:,6)]).^4-1,'linewidth',2,'color','k','linestyle',':');xlim([limit_x_start,limit_x_end])
subplot(2,2,4); plot(0:param.TTT-1,param.bet+dbe_cf,'linewidth',2,'color','k','linestyle',':');xlim([limit_x_start,limit_x_end])

 