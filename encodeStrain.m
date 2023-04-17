function y = encodeStrain( x, resObj )
%encodeStrain encode the chromosome associated to a string x in a
%memory-efficint way. The idea is that the chromosome will roughly have either a
%lot of 1s (when MinCell optimization) or a lot of 0s (usual metabolite
%optimization). I hence remember only the indexes associated to these and
%convert the resulting indexes vector in uint16, that shall give enough
%integer for every metabolic network.
%Input: 
%        -x, i.e. the strain to be encoded
%        -optMode, i.e. the kind of optimization running. Currently
%         supported fields are: i) MinCell. Other to be added yet.
%Output: 
%        -y, i.e. encoded strain.
%Usage:
%        -y = encodeStrain( x, optMode )
%
%**Case MinCell**: since I expect that there will be much more KOs than
%active genes, I remember the active genes only, i.e. y will represents the
%indexes associated to the zeros of x. The remaining entries of x (i.e. the
%ones associated to indexes not appearing in y) are then all ones.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Improvement that can be made to the function: At the cost of adding a    %
%field to y, I could adaptively discern if it is best to encode the       %
%position of the ones, or the position of the zeros. The additional field %
%in y, shall thus indicate which one of the two choices was made. I could %
%also use the minus sign on the indexes to discern that.                  %
%Author: Andrea Patanè                                                    %
%Date: 17/03/2016                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(resObj.flagMinCell )
    y = int16(find(~x));
elseif(resObj.flagKO || resObj.flagRedirector || resObj.flagWC )
     y = int16(find(x));
end
    

end


