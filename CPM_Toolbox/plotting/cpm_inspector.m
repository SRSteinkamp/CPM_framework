function cpm_inspector(pE, pC, M, U, y_true)
% FUNCTION NAME
% cpm_inspector
% 
% DESCRIPTION
% Provides an overview over the underlying model by showing latent and 
% true parameters, and generating data given a set of parameters. 
% Parameters can be interactively changed.
% Correlation and MSE between parameters and true or initial parameters is 
% displayed. Furthermore, the parameters response field is shown.
%
% INPUT:
%   pE - parameter struct of the distribution mean (e.g. prior PRF.M.pE{n} or
%        posterior PRF.Ep{n}
%   pC - struct or covariance matrix of the distribution (e.g. prior
%        PRF.M.pC{n} or posterior PRF.Cp{n}
%   M - the model structure of the (PRF.M)
%   U - the trial information (PRF.U)
%   y_true - real data / data to compare values againts.
% 
% OUTPUT:
%   
% ASSUMPTIONS AND LIMITATIONS:
% Proof of concept.
% 
% TODO:
% * Remove hardcoding of grid parameters
% * Switch between grid parameters
% * allow for 1D grid
% * Fix some visualization issues
% * Clean up code.
%
%
% REVISION HISTORY:
% 29/06/2022 - SRSteinkamp
%   * First implementation. 
%

arguments
   pE
   pC
   M
   U
   y_true = false
end


means = spm_vec(pE);

if strcmpi(class(pC), 'double')

        cov = diag(pC);

elseif isstruct(pC)

        cov = spm_vec(pC);
end

n_params = size(pE, 1);

borders = [means - sqrt(cov) * 3, means,  means + sqrt(cov) * 3]';

if ~y_true
    tmp_y = rwc_simulate_with_prf(spm_unvec(means, pE), M,  U);
else
    tmp_y = y_true;
end
native_vals = cpm_get_true_parameters(spm_unvec(means, pE), M, U);
x     = 1 : 1 : M.ns;


recep  = feval(M.IS, pE, M, U, 'get_response', U(1).gridpoints);
recep  = reshape(recep, U(1).grid.tau(3), U(1).grid.eta(3));
% Open Figure
FigH = figure('units', 'normalized', 'Name', 'Inspector');

% Add parameter selector to change
ParamSelector = uicontrol(FigH,'Style','popupmenu', 'units', 'normalized', ...
    'FontUnits','normalized');
ParamSelector.Position = [0.05 0.35 0.25 0.25];
ParamSelector.String = fieldnames(pE);
ParamSelector.Callback = @selection;
ParamSelector.FontSize = 0.1;

% Slider to change parameter values
SliderH1 = uicontrol(FigH, 'style','slider',...
    'min', borders(1, 1), 'max', borders(3, 1), 'units', 'normalized', ...
    'fontunits','normalized');
SliderH1.Position = [0.4 0.55 0.2 0.05];
SliderH1.FontSize = 0.075;
SliderH1.SliderStep = [0.01 0.1];
set_slider_ticks()

% Text for Parameters
TextParams = uicontrol(FigH, 'style','text',...
     'units', 'normalized', 'fontunits','normalized');
TextParams.Position = [0.05 0.7 0.35 0.30];
TextParams.FontSize = 0.075;
TextParams.HorizontalAlignment = 'left';

TextParams.String = grid_string(ParamSelector.String, means, spm_vec(native_vals));
TextParams.BackgroundColor = [1 1 1];

MainPlot = axes('XLim', [0 M.ns], 'units','normalized', ...
     'NextPlot', 'add', 'position', [0.05, 0.05, 0.9, 0.4]);

ResidPlot = axes('units','normalized', ...
                 'NextPlot', 'add', 'position', [0.725, 0.65, 0.25, 0.3]);


ReceptivePlot = axes('units','normalized', ...
                 'NextPlot', 'add', 'position', [0.45, 0.65, 0.25, 0.3]);
% First plot of tmp_y
LineH = plot(x, [tmp_y, tmp_y], 'parent', MainPlot);
title(MainPlot, "Generated BOLD and prior/provided data")
legend(MainPlot, 'Prior/Provided', 'Generated BOLD')

% Residual Plot
ScatterH = scatter(LineH(1).YData, LineH(2).YData, 'parent', ResidPlot);
title(ResidPlot, ['Generated BOLD vs' newline 'prior/provided data'])

CorrText = text(0.1, 0.9, sprintf('Correlation: %04.3f', corr(LineH(1).YData', LineH(2).YData')), ...
                'Units', 'normalized', 'fontunits', 'normalized', 'fontsize', 0.05, ...
            'parent', ResidPlot);
MSEText = text(0.1, 0.85, sprintf('MSE: %04.3f', mean((LineH(1).YData' - LineH(2).YData').^2)),  ...
                 'Units', 'normalized', 'fontunits', 'normalized', 'fontsize', 0.05, ...
                 'parent', ResidPlot);

% Receptive field
ImageSc= imagesc(U(1).grid.eta(1 : 2), U(1).grid.tau(1 : 2),recep, 'parent', ReceptivePlot);
ReceptivePlot.YLim = U(1).grid.tau(1 : 2);
ReceptivePlot.XLim = U(1).grid.eta(1 : 2);


SliderVal =  uicontrol(FigH, 'style','text',...
   'units', 'normalized', 'fontunits','normalized');
SliderVal.Position = [SliderH1.Position(1) - 0.08, SliderH1.Position(2), 0.075 0.05];
SliderVal.FontSize = 0.75;
SliderVal.BackgroundColor = [1 1 1];
SliderVal.HorizontalAlignment = 'left';

SliderVal.String = num2str(SliderH1.Value);


addlistener(SliderH1, 'Value', 'PostSet', @callbackfn);
addlistener(ParamSelector, 'Value', 'PostSet', @selection);

movegui(FigH, 'center')

function callbackfn(source, eventdata)
    num          = get(eventdata.AffectedObject, 'Value');

    means(ParamSelector.Value) = num;
    tmp_res = rwc_simulate_with_prf(spm_unvec(means, pE), M,  U);
    LineH(2).YData  = tmp_res;

    native_vals = cpm_get_true_parameters(spm_unvec(means, pE), M, U);
    disp(native_vals)
    TextParams.String = grid_string(ParamSelector.String, means, spm_vec(native_vals));
    ScatterH.XData = LineH(1).YData;
    ScatterH.YData = LineH(2).YData;
    SliderVal.String = num2str(round(SliderH1.Value, 3));

    CorrText.String = sprintf('Correlation: %04.3f', corr(LineH(1).YData', LineH(2).YData'));
    MSEText.String = sprintf('MSE: %04.3f', mean((LineH(1).YData' - LineH(2).YData').^2));

    recep  = feval(M.IS, spm_unvec(means, pE), M, U, 'get_response', U(1).gridpoints);
    recep  = reshape(recep, U(1).grid.tau(3), U(1).grid.eta(3));
    ImageSc.CData = recep;
end

function selection(src, event)
    val = ParamSelector.Value;
    str = ParamSelector.String;
    str{val};
    SliderH1.Min = borders(1, ParamSelector.Value);
    SliderH1.Max = borders(3, ParamSelector.Value);
    SliderH1.Value = means(ParamSelector.Value);
    set_slider_ticks()
end


function out_string = grid_string(fields, latent_vals, native_vals)

    max_char = max([length('Parameters'), max(cellfun('length', fields))]);

    out_string = [pad('Parameters', max_char) ' |  latent |  native ' newline];
    for ii = 1 : length(fields)
        tmp_string = sprintf('%s | %+07.4f | %+07.4f', pad(fields{ii}, max_char), ...
                             latent_vals(ii), native_vals(ii));
        out_string = [out_string newline tmp_string];
    end

end

function set_slider_ticks()
   stlblparams = {	'Parent',FigH,...
                    'Units','normalized',...
                    'FontUnits', 'normalized', ...
                    'FontSize', 0.5,...
                    'HorizontalAlignment','center',...
                    'Style','text'};
    ps = SliderH1.InnerPosition;

    lblmargin = 0.05;
    nticks = 7;
    Plbl = [0 ps(2)-0.025 0.05 0.025];
    lblvals = linspace(SliderH1.Min, SliderH1.Max, nticks);
    xpos = linspace(ps(1) + 0.01, ps(1) + ps(3) - 0.01, nticks);
    for k = 1:nticks

        Plbl(1) = xpos(k) - lblmargin /2;

        uicontrol(stlblparams{:},...
            'Position',Plbl,...
            'String',sprintf('%.0f',lblvals(k)),...
            'Tag','slider1ticklabel');
    end

end

  end