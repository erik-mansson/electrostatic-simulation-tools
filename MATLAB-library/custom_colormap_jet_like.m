% Colormap based on 'jet' but customized to have more uniform step contrast (at least for colour-blind).
% Currently the length is fixed at 76 colors.
function cm = custom_colormap_jet_like()

%%
cm = jet(256);

%cm = [cm([1:4:41 44:3:100 102:4:154 157:2:200 202:4:254],:); cm(256,:)*0.95];
cm = [cm([1:4:41 44:3:100 102:4:154 157:2:200 202:4:253],:); cm(256,:)*0.98];
cm([37 41 44 71],:) = [];

% colormap(cm); set(gca,'CLim',[0,length(cm)])