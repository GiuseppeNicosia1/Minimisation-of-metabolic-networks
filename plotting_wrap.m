function plotting_wrap(resObj)
%plotting_wrap, method of the class optResult. It produces plots of the
%result of the optimization/post processing analysis.

labels = computing_labels(resObj);
epsilon = resObj.epsilon;
close all
for ii = 1:(resObj.M/2)
    %%%%%%%%%%%%%%%
    %%MATLAB PLOT%%
    %%%%%%%%%%%%%%%
    figure
    hold on 
    grid on
    if resObj.flagCorner
        title(['Optimization and pareto-analysis results at Corner: ', int2str(ii)]);
    else
        title('Optimization and pareto-analysis results:')
    end
    xlabel(labels{1});
    ylabel(labels{2});
    bioIndex = 2*ii-1;
    synthIndex = 2*ii;
    H = plot( abs(resObj.fit(:,bioIndex)), abs(resObj.fit(:,synthIndex)) ,'b*' );
    H1 = plot( resObj.pfFit(:,bioIndex), resObj.pfFit(:,synthIndex), 'r*' );
    H2 = plot(resObj.bestBioFit(1,ii) , resObj.bestBioFit(2,ii) , 'og' );
    H3 = plot(resObj.bestSynthFit(1,ii) , resObj.bestSynthFit(2,ii), 'ob');
    H4 = plot(resObj.Close2UtopianFit(1,ii) , resObj.Close2UtopianFit(2,ii), 'or');
%     H5 = plot(resObj.bestTradeFit_allCorners(bioIndex) , resObj.bestTradeFit_allCorners(synthIndex), 'xb');
%     H6 = plot(resObj.utopianBio(ii) , resObj.utopianSynth(ii), 'xr'  );
    len = length(epsilon);
    if len > 0
        symb = generate_symbols(len); 	%%genero i simboli per il plot degli epsilon dominate
    end
    K = [];
    for jj=1:len
        aux1 = resObj.epsilonFrontFit{:,1,jj};
        if(~isempty(aux1))
            auxx = aux1(:,bioIndex);
            auxy = aux1(:,synthIndex);
            k = plot( auxx , auxy, symb{jj} );
            K = [K;k];
        else
            epsilon(jj) = -1 ; %flag for empty fronts
        end            
    end
    epsilon(epsilon == -1) = [];   
    legends = generate_legend(epsilon,resObj);
    legend([H;H1;H2;H3;H4;H5;H6;K], legends, 'Location' , 'SouthWest');
    if(resObj.flagCorner)
        savefig([resObj.results_folder, 'postProcessing_Corner_',int2str(ii),'.fig']);
    else
        savefig([resObj.results_folder, 'postProcessing.fig']);
    end
    close all
end

end



function labels=computing_labels(resObj)
	labels{1} = 'biomass [h^{-1}]';
    if(resObj.flagMinCell || resObj.flagWC)
        labels{2} = 'Genes KO';
    elseif(resObj.flagKO || resObj.flagRedirector)
        labels{2} = 'Best Synth';
    end
end


function legends = generate_legend(epsilon,resObj)

	legends{1} = 'Feasible strains';
	legends{2} = 'Pareto strains';
    legends{3} = 'Best Biomass';
    if(resObj.flagMinCell || resObj.flagWC)
        legends{4} = 'Best KOs';
    elseif(resObj.flagKO || resObj.flagRedirector)
        legends{4} = 'Best Synth';
    end
    legends{5} = 'Closest2Utopian';
    legends{6} = 'Closest2Utopian Among All Corners';
    legends{7} = 'Utopian Strain';
	for ii=1:length(epsilon)
		legends{7+ii}=['epsilon strains, epsilon=', num2str(epsilon(ii))];
	end

end

function symb=generate_symbols(len)

	ch={'g','c','m','y','b'};
	ch1={'*','+','o'};

	for ii=1:len 
		symb{ii}=[ch{mod(ii,5)+1} , ch1{mod(ii,3)+1} ];
	end






end

