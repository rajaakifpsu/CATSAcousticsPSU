%Purpose: to create a power spectral density plot for a single mic

%Inputs: fpwelch is vector or frequencies corresponding to PSD
%       PSD is the matrix representing power spectral density, where each
%       row is a frequency, each column is a mic.
%       graphtitle is the title of the graph it makes
%       micNum is/are the number(s) of the mics to use. If multiple
%       numbers in a vector, then those mics will all be plotted on the same plot.
%       xlimits is a 2-element increasing vector of the limits of the x
%       axis in the plot.
%       ylimits the same but for y limits.
%       graphPos is an array of 4 integers indicating the x and y location
%       and x and y size of the graph to be made.
%       newfig is boolean for if a new figure is to be created. 

%Outputs: a figure

function createPowerSpectralDensityPlotSingular(fpwelch, PSD, graphtitle, micNum, xlimits, ylimits, graphPos, newfig)
    if newfig %only creates new figure if newfig true
        figure()
%     else
%         hold on;
    end
    plot(fpwelch,PSD(:,micNum))
    xlabel('Frequency (Hz)')
    ylabel('Power Spectral Density (dB/Hz)')
    xlim(xlimits)
    ylim(ylimits)
    title(graphtitle)
    set(gcf,'position',graphPos)
    grid on
end