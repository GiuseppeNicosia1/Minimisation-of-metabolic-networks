function x = decodeStrain( y , resObj )
%decodeStrain decode a strain from the compressed form to the usual form 
%(namely that one used during the optimization).
%Input:
%        -y, i.e. the encoded strain
%        -option, i.e. the usual structure used thorought the optimization
%Output:
%        -x, i.e. the decode strain
%Usage:
%        -x = decodeStrain( y , option );
%This function is the "inverse" of encodeStrain. In the sense that both the
%following statements are true:
% isequal(x,  decodeStrain( encodeStrain(x,option.optMode) , option) )
% isequal(y , encodeStrain( decodeStrain(y,option) , option.optMode) )
%in which I use the same structure "option" to assure that I am of course
%talking about the same optimization.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Improvements that may be made to the function: see encodeStrain.         %
%Author: Andrea Patanè                                                    %
%Date: 18/03/2016                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(resObj.flagMinCell)
    x = ones(1,resObj.V);
    x(y) = 0;
elseif(resObj.flagKO || resObj.flagRedirector || resObj.flagWC)
    x = zeros(1,resObj.V);
    x(y) = 1;
end

end


