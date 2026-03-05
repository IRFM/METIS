function hf = metis4muscle3_layout(internal_state)


% save contexte_layout

% create a figure
hf = figure('tag','metis4muscle3_actor_layout','name','METIS actor for MUSCLE3','defaultaxesfontsize',12,'defaultaxesfontweight','bold', ...
    'defaultaxesfontname','times','defaultlinelinewidth',1,'color',[1 1 1], 'toolbar','figure','position',[100,100,1440,960]);

% actor rectangle
layout.x0 = -15;
layout.y0 = 0;
layout.dx = 150;
layout.dy = 100;
layout.axes_margin_left  = -5;
layout.axes_margin_right = 30;
layout.axes_margin_up = 30;
layout.axes_margin_down = -45;
layout.offset_x = 4;
layout.offset_y = 8;

rectangle('Position',[layout.x0, layout.y0, layout.dx, layout.dy],'Curvature',0.05);
hold on
% make axes invisible
set(gca,'xtick',[]);
set(gca,'xcolor',[1 1 1])
set(gca,'ytick',[]);
set(gca,'ycolor',[1 1 1])
set(gca,'xlim',[layout.x0 + layout.axes_margin_left ,layout.x0 + layout.dx + layout.axes_margin_right]);
set(gca,'ylim',[layout.y0 + layout.axes_margin_down,layout.y0 + layout.dy + layout.axes_margin_up]);
ha_main = gca;

% position pointers
layout.init.y = layout.y0 + layout.dy + 1;
layout.init.dx = 4;
layout.init.x = layout.offset_x;

layout.input_f.x= layout.x0-1;
layout.input_f.dy = 4;
layout.input_f.y= layout.y0 + layout.dy - 3 .* layout.offset_y;

layout.input_s.x= layout.x0 -1;
layout.input_s.dy = 4;
layout.input_s.y= layout.y0 + layout.dy / 2 - layout.offset_y;

layout.output_f.x = layout.x0 + layout.dx + 1;
layout.output_f.dy = 3.33;
layout.output_f.y = layout.y0 + layout.dy - layout.offset_y / 2;

layout.output_i.x = layout.x0 + layout.dx + 1;
layout.output_i.dy = 3.33;
layout.output_i.y = layout.y0 + layout.dy / 2 - layout.offset_y / 2;

layout.setting.y = layout.y0 - 1;
layout.setting.dx = 4;
layout.setting.x = layout.offset_x;

layout.option.y = layout.y0 + layout.dy + 1;
layout.option.dx = 4;
layout.option.x = layout.x0 + layout.dx -2 * layout.offset_x;

layout.args.y  = layout.y0 - 1;
layout.args.dx = 4;
layout.args.x = layout.x0 + layout.dx -2 * layout.offset_x;

layout.marker_size = 6;
layout.text_margin = 3;

% to compute barycenter
layout_mem = layout;

% loop on connected ports
noms = fieldnames(internal_state.actor.io.status);
for k=1:length(noms)
    if internal_state.actor.io.status.(noms{k})
        ck = [0 0 0];
    else
        ck = [0.7 0.7 0.7];
    end
    % switch on direction
    switch internal_state.actor.io.list.(noms{k}).direction
        case 'output'
            switch internal_state.actor.io.list.(noms{k}).operator
                case 'O_F'
                    plot(layout.output_f.x,layout.output_f.y,'color','k','linestyle','none','marker','d','markersize',layout.marker_size,'MarkerEdgeColor',ck,'MarkerFaceColor',ck);
                    text(layout.output_f.x + layout.text_margin,layout.output_f.y,noms{k},'color',ck,'Interpreter','none','Rotation',0,'HorizontalAlignment','left');
                    layout.output_f.y = layout.output_f.y - layout.output_f.dy;
                case 'O_I'
                    plot(layout.output_i.x,layout.output_i.y,'color','k','linestyle','none','marker','o','markersize',layout.marker_size,'MarkerEdgeColor',ck,'MarkerFaceColor',ck);
                    text(layout.output_i.x + layout.text_margin,layout.output_i.y,noms{k},'color',ck,'Interpreter','none','Rotation',0,'HorizontalAlignment','left');
                    layout.output_i.y = layout.output_i.y - layout.output_i.dy;
                otherwise
                    error('unknown operator');
            end
        case 'input'
            switch internal_state.actor.io.list.(noms{k}).operator
                case 'F_INIT'
                    plot(layout.input_f.x,layout.input_f.y,'color',ck,'linestyle','none','marker','d','markersize',layout.marker_size,'MarkerEdgeColor',ck,'MarkerFaceColor','none');
                    text(layout.input_f.x-layout.text_margin,layout.input_f.y,noms{k},'color',ck,'Interpreter','none','Rotation',0,'HorizontalAlignment','right');
                    layout.input_f.y = layout.input_f.y - layout.input_f.dy;
                case 'S'
                    plot(layout.input_s.x,layout.input_s.y,'color',ck,'linestyle','none','marker','o','markersize',layout.marker_size,'MarkerEdgeColor',ck,'MarkerFaceColor','none');
                    text(layout.input_s.x-layout.text_margin,layout.input_s.y,noms{k},'color',ck,'Interpreter','none','Rotation',0,'HorizontalAlignment','right');
                    layout.input_s.y = layout.input_s.y - layout.input_s.dy;
                otherwise
                    error('unknown operator');
            end
        case 'init'
            switch internal_state.actor.io.list.(noms{k}).operator
                case 'F_INIT'
                    plot(layout.init.x,layout.init.y,'color',ck,'linestyle','none','marker','d','markersize',layout.marker_size,'MarkerEdgeColor',ck,'MarkerFaceColor','none');
                    text(layout.init.x,layout.init.y+layout.text_margin,noms{k},'color',ck,'Interpreter','none','Rotation',-90,'HorizontalAlignment','right');
                    layout.init.x = layout.init.x + layout.init.dx;
                case 'S'
                    error('not yet implemented');
                otherwise
                    error('unknown operator');
            end
            
        otherwise
            error('unknown direction');
    end
end

% legend for I/O
text((layout.init.x + layout_mem.init.x) / 2 ,layout.init.y - 1.5 * layout.text_margin,'Initialization ports (F_INIT)', ...
      'color','k','Interpreter','none','Rotation',0,'HorizontalAlignment','center','VerticalAlignment','middle','EdgeColor','k');
text(layout.input_f.x + 1.5 * layout.text_margin,(layout.input_f.y + layout_mem.input_f.y) / 2,'Input ports (F_INIT)', ...
      'color','k','Interpreter','none','Rotation',-90,'HorizontalAlignment','center','VerticalAlignment','middle','EdgeColor','k');
text(layout.input_s.x + 1.5 * layout.text_margin,(layout.input_s.y + layout_mem.input_s.y) / 2,'Input ports (S)', ...
      'color','k','Interpreter','none','Rotation',-90,'HorizontalAlignment','center','VerticalAlignment','middle','EdgeColor','k');
text(layout.output_f.x - 1.5 * layout.text_margin,(layout.output_f.y + layout_mem.output_f.y) / 2,'Output ports (O_F)', ...
      'color','k','Interpreter','none','Rotation',-90,'HorizontalAlignment','center','VerticalAlignment','middle','EdgeColor','k');
text(layout.output_i.x - 1.5 * layout.text_margin,(layout.output_i.y + layout_mem.output_i.y) / 2,'Output ports (O_I)', ...
      'color','k','Interpreter','none','Rotation',-90,'HorizontalAlignment','center','VerticalAlignment','middle','EdgeColor','k');
  
% setting
for k=1:length(internal_state.actor.settings.list)
  if ~isempty(strmatch(internal_state.actor.settings.list{k}{1},internal_state.actor.settings.value.instance_listofsettings,'exact'))
      val_def = internal_state.actor.settings.list{k}{3};
      val     = eval(sprintf('internal_state.actor.settings.value.%s',strrep(internal_state.actor.settings.list{k}{1},'_','.')));
      if ischar(val)
          if strcmp(val,val_def)
              ck = [0 0.447 0.741];
          elseif isempty(val)
              ck = [0.7 0.7 0.7];              
          else
              ck = [0 0 0];
          end
      elseif val_def == val
          ck = [0 0.447 0.741];
      elseif val == 0
          ck = [0.929 0.6940 0.125];
      else
          ck = [0.466 0.674 0.188];          
      end
  else
          ck = [0 0.447 0.741];
  end
  plot(layout.setting.x,layout.setting.y,'color',ck,'linestyle','none','marker','*','markersize',layout.marker_size,'MarkerEdgeColor',ck,'MarkerFaceColor','none');
  text(layout.setting.x,layout.setting.y  - layout.text_margin,internal_state.actor.settings.list{k}{1},'color',ck,'Interpreter','none','Rotation',-90,'HorizontalAlignment','left');
  layout.setting.x = layout.setting.x + layout.setting.dx;
end
text((layout.setting.x + layout_mem.setting.x) / 2 ,layout.setting.y + 1.5 * layout.text_margin,'Actor settings', ...
      'color','k','Interpreter','none','Rotation',0,'HorizontalAlignment','center','VerticalAlignment','middle','EdgeColor','k');


% METIS options
if isfield(internal_state.actor.settings.value,'metis_option')
  noms = fieldnames(internal_state.actor.settings.value.metis_option);
  option = metis4imas(1);
  for k=1:length(noms)
      val = internal_state.actor.settings.value.metis_option.(noms{k});
      % for convinience and to be user friendly, metis parameter
      metis_name = noms{k};
      % can be design by is label
      if ~isfield(option.valeur,metis_name)
          % search for label
          nomsl = fieldnames(option.label);
          for l=1:length(nomsl)
              switch  option.label.(nomsl{l})
                  case metis_name
                      metis_name = nomsl{l};
                      break
              end
          end
      end
      if isfield(option.valeur,metis_name)
          val_def = option.valeur.(metis_name);
          if ischar(val)
              if strcmp(val,val_def)
                  ck = [0 0.447 0.741];
              elseif isempty(val)
                  ck = [0.7 0.7 0.7];
              else
                  ck = [0 0 0];
              end
          elseif val_def == val
              ck = [0 0.447 0.741];
          elseif val == 0
              ck = [0.929 0.6940 0.125];
          else
              ck = [0.466 0.674 0.188];
          end
      else
          % error
          ck = [0.850 0.325 0.098];
      end
      plot(layout.option.x,layout.option.y,'color',ck,'linestyle','none','marker','h','markersize',layout.marker_size,'MarkerEdgeColor',ck,'MarkerFaceColor','none');
      text(layout.option.x,layout.option.y+layout.text_margin,noms{k},'color',ck,'Interpreter','none','Rotation',-90,'HorizontalAlignment','right');
      layout.option.x = layout.option.x - layout.option.dx;
  end
  text(layout_mem.option.x - 4 * layout.option.dx - layout.text_margin,layout.option.y - 1.5 * layout.text_margin,'METIS parameters (overwritten)', ...
      'color','k','Interpreter','none','Rotation',0,'HorizontalAlignment','center','VerticalAlignment','middle','EdgeColor','k');
else
  %text((layout.option.x + layout_mem.option.x) / 2 - layout.text_margin,layout.option.y - 1.5 * layout.text_margin,'METIS parameters (None defined)', ...
  %    'color','k','Interpreter','none','Rotation',0,'HorizontalAlignment','center','VerticalAlignment','middle','EdgeColor','k');
  text(layout_mem.option.x - 4 * layout.option.dx - layout.text_margin,layout.option.y - 1.5 * layout.text_margin,'METIS parameters (None defined)', ...
      'color','k','Interpreter','none','Rotation',0,'HorizontalAlignment','center','VerticalAlignment','middle','EdgeColor','k');
end

  
filename = internal_state.actor.settings.value.metis.filename;
[~,filename] = fileparts(filename);
filename = fullfile('...',filename);
text(layout.x0 + 0.5.* layout.dx , layout.option.y - 1.5 * layout.text_margin,sprintf('METIS file: %s',filename), ...
       'color','k','Interpreter','none','Rotation',0,'HorizontalAlignment','center','VerticalAlignment','middle','EdgeColor','k','fontsize',8);
   
text(layout.x0 + 0.5.* layout.dx , layout.option.y + 5 * layout.text_margin,'METIS actor for MUSCLE3', ...
       'color','k','Interpreter','none','Rotation',0,'HorizontalAlignment','center','VerticalAlignment','middle', ...
       'EdgeColor','k','fontsize',14,'FontWeight','bold');
  

% Actor arguments
if internal_state.workflow.checkpoint_on_off
    ck = [0.466 0.674 0.188];
else
    ck = [0.850 0.325 0.098];
end
plot(layout.args.x,layout.args.y,'color',ck,'linestyle','none','marker','pentagram','markersize',layout.marker_size,'MarkerEdgeColor',ck,'MarkerFaceColor','none');
text(layout.args.x,layout.args.y  - layout.text_margin,'checkpoint_on_off','color',ck,'Interpreter','none','Rotation',-90,'HorizontalAlignment','left');
layout.args.x = layout.args.x - layout.args.dx;
if isfield(internal_state.workflow,'overlay_on_off')  && internal_state.workflow.overlay_on_off
    ck = [0.466 0.674 0.188];
else
    ck = [0.850 0.325 0.098];
end
plot(layout.args.x,layout.args.y,'color',ck,'linestyle','none','marker','pentagram','markersize',layout.marker_size,'MarkerEdgeColor',ck,'MarkerFaceColor','none');
text(layout.args.x,layout.args.y  - layout.text_margin,'overlay_on_off','color',ck,'Interpreter','none','Rotation',-90,'HorizontalAlignment','left');
layout.args.x = layout.args.x - layout.args.dx;
if internal_state.workflow.profiler_status
    ck = [0.466 0.674 0.188];
else
    ck = [0.850 0.325 0.098];
end
plot(layout.args.x,layout.args.y,'color',ck,'linestyle','none','marker','pentagram','markersize',layout.marker_size,'MarkerEdgeColor',ck,'MarkerFaceColor','none');
text(layout.args.x,layout.args.y  - layout.text_margin,'Matlab profiler','color',ck,'Interpreter','none','Rotation',-90,'HorizontalAlignment','left');
layout.args.x = layout.args.x - layout.args.dx;
switch internal_state.workflow.debug_level
    case 0
        ck = [0.850 0.325 0.098];       
    case 1
        ck = [0.466 0.674 0.188];
    otherwise
        ck = [0 0.447 0.741];
end
plot(layout.args.x,layout.args.y,'color',ck,'linestyle','none','marker','pentagram','markersize',layout.marker_size,'MarkerEdgeColor',ck,'MarkerFaceColor','none');
text(layout.args.x,layout.args.y  - layout.text_margin,'Debug level','color',ck,'Interpreter','none','Rotation',-90,'HorizontalAlignment','left');
layout.args.x = layout.args.x - layout.args.dx;
text((layout.args.x + layout_mem.args.x) / 2 ,layout.args.y + 1.5 * layout.text_margin,'Actor arguments', ...
      'color','k','Interpreter','none','Rotation',0,'HorizontalAlignment','center','VerticalAlignment','middle','EdgeColor','k');
    


% clock
switch internal_state.actor.settings.value.metis.clock.internal
  case 'on'
      % display internal clock
      clock_mode = 'internal';
  otherwise
      clock_mode = 'external';
end
pos = get(ha_main,'position');
hck = axes('position',[pos(1) + 0.6 .* pos(3) ,pos(2) + 0.58 .* pos(4) , pos(3) / 6, 0.15 .* pos(4)],'Box','on', ...
    'XTick','','XTickLabel','','XTickLabelMode','manual','YTick','','YTickLabel','','YTickLabelMode','manual', ...
    'Xlim',[0 1],'Ylim',[0,5]);
title('Clock','EdgeColor','k');
text(0.05,1.1 .* (4 -0.1) ,sprintf('mode: %s',clock_mode), 'color','k','Interpreter','none','Rotation',0,'HorizontalAlignment','left','VerticalAlignment','middle','EdgeColor','none','fontsize',8);
text(0.05,1.1 .* (3 -0.1) ,sprintf('start: %6.3g (s)',internal_state.actor.settings.value.metis.clock.time.start), 'color','k','Interpreter','none','Rotation',0,'HorizontalAlignment','left','VerticalAlignment','middle','EdgeColor','none','fontsize',8);
text(0.05,1.1 .* (2 -0.1) ,sprintf('stop: %6.3g (s)',internal_state.actor.settings.value.metis.clock.time.stop), 'color','k','Interpreter','none','Rotation',0,'HorizontalAlignment','left','VerticalAlignment','middle','EdgeColor','none','fontsize',8);
text(0.05,1.1 .* (1 -0.1) ,sprintf('time step: %6.3g (s)',internal_state.actor.settings.value.metis.clock.time.dt), 'color','k','Interpreter','none','Rotation',0,'HorizontalAlignment','left','VerticalAlignment','middle','EdgeColor','none','fontsize',8);

% state: display non zero state
pos = get(ha_main,'position');
noms = fieldnames(internal_state.workflow.state);
hstate = axes('position',[pos(1) + 0.6 .* pos(3) ,pos(2) + 0.33 .* pos(4) , pos(3) / 6, 0.2 .* pos(4)],'Box','on', ...
    'XTick','','XTickLabel','','XTickLabelMode','manual','YTick','','YTickLabel','','YTickLabelMode','manual', ...
    'Xlim',[0 1],'Ylim',[0,length(noms) + 1]);
title('State','EdgeColor','k');
for k=1:length(noms)
    if ~ischar(internal_state.workflow.state.(noms{k}))
        
        if internal_state.workflow.state.(noms{k})
            ck = [0.466 0.674 0.188];
        else
            ck = [0.850 0.325 0.098];
        end
        text(0.05,1.1 .* (k -0.1) ,noms{k}, 'color',ck,'Interpreter','none','Rotation',0,'HorizontalAlignment','left','VerticalAlignment','middle','EdgeColor','none','fontsize',8);
    end
end



% external data
%search for connected ports
noms = fieldnames(internal_state.actor.io.status);
ext_equi = false;
ext_cp = false;
ext_cs = false;
for k=1:length(noms)
    if internal_state.actor.io.status.(noms{k})
       if (startsWith(noms{k},'plasma_profiles') || startsWith(noms{k},'core_profiles'))  && ...
          (isempty(strfind(noms{k},'_output_')) || isempty(strfind(noms{k},'_out_'))) && ...
           internal_state.actor.io.status.(noms{k})    
           ext_cp = true;
       end
       if (startsWith(noms{k},'plasma_sources') || startsWith(noms{k},'core_sources'))  && ...
          (isempty(strfind(noms{k},'_output_')) || isempty(strfind(noms{k},'_out_'))) && ...
           internal_state.actor.io.status.(noms{k})
           ext_cs = true;
       end
       if startsWith(noms{k},'equilibrium') && ...
          (isempty(strfind(noms{k},'_output_')) || isempty(strfind(noms{k},'_out_'))) && ...
          internal_state.actor.io.status.(noms{k})
           ext_equi = true;
       end
    end
end

pos = get(ha_main,'position');
noms = fieldnames(internal_state.actor.settings.value.metis_external_data);
hex = axes('position',[pos(1) + 0.07 .* pos(3) ,pos(2) + 0.32 .* pos(4) , 0.24 * pos(3), 0.4 .* pos(4)],'Box','on', ...
           'XTick','','XTickLabel','','XTickLabelMode','manual','YTick','','YTickLabel','','YTickLabelMode','manual', ...
           'Xlim',[0 1],'Ylim',[0,length(noms) + 2]);
title('External data','EdgeColor','k');


for k=1:length(noms)
    leg_loc = '';
    if ext_equi
       ck = [0.466 0.674 0.188];
       leg_loc = 'with equi';
    else
       ck = [0 0.447 0.741];        
       leg_loc = 'w/o equi';
    end
    if ~internal_state.actor.settings.value.metis_external_data.(noms{k}).state
        ck = [0.7 0.7 0.7];
        leg_loc = 'offline';
    else
        switch noms{k}
            case  {'equilibrium','current_diffusion'}
                if ~ext_equi
                   ck = [0.850 0.325 0.098]; 
                   leg_loc = 'equilibrium missing';
               end
            case {'electron_temperature','electron_density','ion_temperature','charge_effective','data_momentum'}
                if ~ext_cp
                   ck = [0.850 0.325 0.098]; 
                   leg_loc = 'plasma_profiles missing';
               end
            case {'LHCD','ECRH','ICRH','NBI','PLINE','runaway'}
                if ~ext_cs
                   ck = [0.850 0.325 0.098]; 
                   leg_loc = 'plasma_sources missing';
                end
        end
    end        
    text(0.05,1.1 .* (length(noms) - k + 1 -0.1) ,sprintf('%s (%s)',noms{k},leg_loc), 'color',ck,'Interpreter','none','Rotation',0,'HorizontalAlignment','left','VerticalAlignment','middle','EdgeColor','none','fontsize',8);
end



% overlay METIS image
pos = get(ha_main,'position');
him = axes('position',[pos(1) + pos(3) / 3,pos(2) + pos(4) / 3.3, pos(3) / 4, pos(4) / 2.2]);
% read image
im = imread(fullfile(fileparts(which('metis')),'doc','METIS_layout.png'));
image(im)
set(him,'xtick',[]);
set(him,'xcolor',[1 1 1])
set(him,'ytick',[]);
set(him,'ycolor',[1 1 1])


% add legends
% for ports
pos = get(ha_main,'position');
hl1 = axes('position',[pos(1) - 0.15 * pos(3) ,pos(2) +  0.87 * pos(4) , pos(3) / 6, 0.1 .* pos(4)],'Box','on', ...
    'XTick','','XTickLabel','','XTickLabelMode','manual','YTick','','YTickLabel','','YTickLabelMode','manual', ...
    'Xlim',[0 1],'Ylim',[0,5]);
title('Ports','EdgeColor','k','FontWeight','normal','fontsize',10);
text(0.05,1.1 .* (4 -0.1) ,'* connected', 'color','k','Interpreter','none','Rotation',0,'HorizontalAlignment','left','VerticalAlignment','middle','EdgeColor','none','fontsize',8);
text(0.05,1.1 .* (3 -0.1) ,'* disconnected', 'color',[0.7 0.7 0.7],'Interpreter','none','Rotation',0,'HorizontalAlignment','left','VerticalAlignment','middle','EdgeColor','none','fontsize',8);
text(0.05,1.1 .* (1 -0.1) ,sprintf('S port data overwite\nF_INIT port data'), 'color','k','Interpreter','none','Rotation',0,'HorizontalAlignment','left','VerticalAlignment','middle','EdgeColor','none','fontsize',8);


% for setting
pos = get(ha_main,'position');
hl2 = axes('position',[pos(1) - 0.15 * pos(3) ,pos(2) +  0 * pos(4) , pos(3) / 6, 0.15 .* pos(4)],'Box','on', ...
    'XTick','','XTickLabel','','XTickLabelMode','manual','YTick','','YTickLabel','','YTickLabelMode','manual', ...
    'Xlim',[0 1],'Ylim',[0,6]);
title('Settings','EdgeColor','k','FontWeight','normal','fontsize',10);
text(0.05,1.1 .* (5 -0.1) ,'* empty', 'color',[0.7 0.7 0.7],'Interpreter','none','Rotation',0,'HorizontalAlignment','left','VerticalAlignment','middle','EdgeColor','none','fontsize',8);
text(0.05,1.1 .* (4 -0.1) ,'* default', 'color',[0 0.447 0.741],'Interpreter','none','Rotation',0,'HorizontalAlignment','left','VerticalAlignment','middle','EdgeColor','none','fontsize',8);
text(0.05,1.1 .* (3 -0.1) ,'* customized (string)', 'color',[0 0 0],'Interpreter','none','Rotation',0,'HorizontalAlignment','left','VerticalAlignment','middle','EdgeColor','none','fontsize',8);
text(0.05,1.1 .* (2 -0.1) ,'* true or customized (numeric)', 'color',[0.466 0.674 0.188],'Interpreter','none','Rotation',0,'HorizontalAlignment','left','VerticalAlignment','middle','EdgeColor','none','fontsize',8);
text(0.05,1.1 .* (1 -0.1) ,'* false or zero (numeric)', 'color',[0.929 0.6940 0.125],'Interpreter','none','Rotation',0,'HorizontalAlignment','left','VerticalAlignment','middle','EdgeColor','none','fontsize',8);

 
% for parameters
pos = get(ha_main,'position');
hl3 = axes('position',[pos(1) + 0.87 * pos(3) ,pos(2) +  0.87 * pos(4) , pos(3) / 6, 0.15 .* pos(4)],'Box','on', ...
    'XTick','','XTickLabel','','XTickLabelMode','manual','YTick','','YTickLabel','','YTickLabelMode','manual', ...
    'Xlim',[0 1],'Ylim',[0,7]);
title('Parameters','EdgeColor','k','FontWeight','normal','fontsize',10);
text(0.05,1.1 .* (6 -0.1) ,'* empty', 'color',[0.7 0.7 0.7],'Interpreter','none','Rotation',0,'HorizontalAlignment','left','VerticalAlignment','middle','EdgeColor','none','fontsize',8);
text(0.05,1.1 .* (5 -0.1) ,'* default', 'color',[0 0.447 0.741],'Interpreter','none','Rotation',0,'HorizontalAlignment','left','VerticalAlignment','middle','EdgeColor','none','fontsize',8);
text(0.05,1.1 .* (4 -0.1) ,'* customized (string)', 'color',[0 0 0],'Interpreter','none','Rotation',0,'HorizontalAlignment','left','VerticalAlignment','middle','EdgeColor','none','fontsize',8);
text(0.05,1.1 .* (3 -0.1) ,'* true or customized (numeric)', 'color',[0.466 0.674 0.188],'Interpreter','none','Rotation',0,'HorizontalAlignment','left','VerticalAlignment','middle','EdgeColor','none','fontsize',8);
text(0.05,1.1 .* (2 -0.1) ,'* false or zero (numeric)', 'color',[0.929 0.6940 0.125],'Interpreter','none','Rotation',0,'HorizontalAlignment','left','VerticalAlignment','middle','EdgeColor','none','fontsize',8);
text(0.05,1.1 .* (1 -0.1) ,'* error', 'color',[0.850 0.325 0.098],'Interpreter','none','Rotation',0,'HorizontalAlignment','left','VerticalAlignment','middle','EdgeColor','none','fontsize',8);
     
% for arguments
pos = get(ha_main,'position');
hl3 = axes('position',[pos(1) + 0.87 * pos(3) ,pos(2) +  0.05 * pos(4) , pos(3) / 6, 0.1 .* pos(4)],'Box','on', ...
    'XTick','','XTickLabel','','XTickLabelMode','manual','YTick','','YTickLabel','','YTickLabelMode','manual', ...
    'Xlim',[0 1],'Ylim',[0,4]);
title('Arguments','EdgeColor','k','FontWeight','normal','fontsize',10);
text(0.05,1.1 .* (3 -0.1) ,'* false ', 'color',[0.850 0.325 0.098],'Interpreter','none','Rotation',0,'HorizontalAlignment','left','VerticalAlignment','middle','EdgeColor','none','fontsize',8);
text(0.05,1.1 .* (2 -0.1) ,'* true', 'color',[0.466 0.674 0.188],'Interpreter','none','Rotation',0,'HorizontalAlignment','left','VerticalAlignment','middle','EdgeColor','none','fontsize',8);
text(0.05,1.1 .* (1 -0.1) ,'* numeric > 1', 'color',[0 0.447 0.741],'Interpreter','none','Rotation',0,'HorizontalAlignment','left','VerticalAlignment','middle','EdgeColor','none','fontsize',8);


% print figure
if nargout == 0
    %    print('METIS_actor_for_MUSCLE3_layout','-dpng');
    try
        % print trigger error on some system
         exportgraphics(hf, 'METIS_actor_for_MUSCLE3_layout.png', 'BackgroundColor', [1 1 1], 'Resolution', 300);   
    catch
        % just save Matlab object that can be reopen for debug in a Matlab
        % session
        savefig(hf,'METIS_actor_for_MUSCLE3_layout');
    end
end



     
