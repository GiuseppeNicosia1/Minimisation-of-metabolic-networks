function minSolsCompGenes(model)

fbamodel = load(['solutions' model filesep 'fbamodel.mat']);
fbamodel = fbamodel.fbamodel;
load(['solutions' model filesep 'minimalSolutions.mat']);

fbamodel.genes = regexprep(fbamodel.genes,'_deleted','');
% fbamodel.genes = regexprep(fbamodel.genes,'_','-');

[M,N] = size(minSol);

KOs = sum(minSol(:,1:fbamodel.ngenes),1);
int = find(KOs == 0);
uni = find(KOs > 0 & KOs < M);
% uni = setdiff(uni,int); % only the difference?
KOs(uni) = M - KOs(uni);
[~, I] = sortrows(KOs(uni)',-1); 
nu = find(KOs == M);

tab = cell(length(int)+1,1);
tab(1,:) = {'Genes Intersection'};
for i = 1:length(int)
%     N = strcmp(fbamodel.genes(int(i)),geneNames(:,2));
    tab(i+1,:) = fbamodel.genes(int(i));
end

tab2 = cell(length(uni)+1,2);
tab2(1,:) = {'Frequency' 'Active Genes Names'};
for i = 1:length(uni)
%     N = strcmp(fbamodel.genes(uni(I(i))),geneNames(:,2));
%     tab2(i+1,1:2+2*any(N)) = {num2str(KOs(uni(I(i)))/M) fbamodel.genes{uni(I(i))} geneNames{N,4} geneNames{N,5}};
%     tab2(i+1,:) = {num2str(KOs(uni(I(i)))/M) fbamodel.genes{uni(I(i))} geneNames{N,4} geneNames{N,5}};
    tab2(i+1,:) = {num2str(KOs(uni(I(i)))/M) fbamodel.genes{uni(I(i))}};
end

tab3 = cell(length(nu)+1,1);
tab3(1,:) = {'Genes Never Active'};
for i = 1:length(nu)
%     N = strcmp(fbamodel.genes(nu(i)),geneNames(:,2));
    tab3(i+1,:) = fbamodel.genes(nu(i));
end

if fbamodel.flagEssGenes
    temp = fbamodel.genesBU(fbamodel.essentialGenes);
    tab4 = cell(length(temp)+1,1);
    tab4(1,:) = {'Essential Genes'};
    for i = 1:length(temp)
%         N = strcmp(temp(i),geneNames(:,2));
        tab4(i+1,:) = temp(i);
    end
end

if fbamodel.flagEssGenes
    save(['solutions' model filesep 'GenesAnTabs.mat'],'tab','tab2','tab3','tab4')
else
    save(['solutions' model filesep 'GenesAnTabs.mat'],'tab','tab2','tab3')
end

% % name = ['excel' filesep 'genesAnalysis' analysis '_' folder(1:end-1) '.xls'];
% % if exist(name,'file')
% %     delete(name)
% % end
% % tab(2:end,:) = sortrows(tab(2:end,:),1);
% % % tab2(2:end,:) = sortrows(tab2(2:end,:),2);
% % tab3(2:end,:) = sortrows(tab3(2:end,:),1);
% % tab4(2:end,:) = sortrows(tab4(2:end,:),1);
% % xlswrite(name,tab,'Active Genes Intersection');
% % xlswrite(name,tab2,'Active Genes Union (noInt)');
% % xlswrite(name,tab3,'Genes always KO');
% % xlswrite(name,tab4,'Essential Genes');
% % newExcel = actxserver('excel.application');
% % excelWB = newExcel.Workbooks.Open([pwd filesep name],0,false);
% % newExcel.Sheets.Item(1).Delete;
% % for i = 1:4
% %     newExcel.Sheets.Item(i).Columns.AutoFit;
% %     set(newExcel.Sheets.Item(i).Rows.Item(1).Font, 'Bold', true);
% % end
% % excelWB.Save();
% % excelWB.Close();
% % newExcel.Quit();
% % delete(newExcel);


end