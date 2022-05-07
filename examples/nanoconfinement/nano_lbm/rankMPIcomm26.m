function [commF,sendRank]=rankMPIcomm26(threadi,threadj,threadk,periodicX_,periodicY_,periodicZ_,numPx,numPy,numPz)
%returns if commumication between ranks exist 'commF' and what the
%communicating rank is 'sendRank' for (26 comm directions), if no comm
%returns -1
                                                                               %corner passes
cx  = [1, -1,  0,  0,  0,  0,  1, -1,  1, -1,  1, -1,  1, -1,  0,  0,  0,  0,  1, -1,  1, -1,  1, -1,  1, -1];
cy  = [0,  0,  1, -1,  0,  0,  1,  1, -1, -1,  0,  0,  0,  0,  1, -1,  1, -1,  1,  1, -1, -1,  1,  1, -1, -1];
cz  = [0,  0,  0,  0,  1, -1,  0,  0,  0,  0,  1,  1, -1, -1,  1,  1, -1, -1   1,  1,  1,  1, -1, -1, -1, -1];

sendRank=-ones(1,26);
if periodicZ_==0 && periodicY_==0 && periodicX_==0
    for pop = 1:26
        threadiN=threadi+cx(pop);
        threadjN=threadj+cy(pop);
        threadkN=threadk+cz(pop);
        commF(pop)=1;
        if threadiN>numPx-1
            commF(pop)=0;
        elseif threadiN<0
            commF(pop)=0;
        end
        if threadjN>numPy-1
            commF(pop)=0;
        elseif threadjN<0
            commF(pop)=0;
        end
        if threadkN>numPz-1
            commF(pop)=0;
        elseif threadkN<0
            commF(pop)=0;
        end
        
        if commF(pop)==1
            sendRank(pop)=threadiN+threadjN*(numPx)+threadkN*(numPy*numPx);
        end
    end
end
if periodicZ_==1 && periodicY_==0 && periodicX_==0
    for pop = 1:26
        threadiN=threadi+cx(pop);
        threadjN=threadj+cy(pop);
        threadkN=threadk+cz(pop);
        commF(pop)=1;
        if threadiN>numPx-1
            commF(pop)=0;
        elseif threadiN<0
            commF(pop)=0;
        end
        if threadjN>numPy-1
            commF(pop)=0;
        elseif threadjN<0
            commF(pop)=0;
        end
        if threadkN>numPz-1
            threadkN=0;
        elseif threadkN<0
            threadkN=numPz-1;
        end
        
        if commF(pop)==1
            sendRank(pop)=threadiN+threadjN*(numPx)+threadkN*(numPy*numPx);
        end
    end
end
if periodicZ_==0 && periodicY_==1 && periodicX_==0
    for pop = 1:26
        threadiN=threadi+cx(pop);
        threadjN=threadj+cy(pop);
        threadkN=threadk+cz(pop);
        commF(pop)=1;
        if threadiN>numPx-1
            commF(pop)=0;
        elseif threadiN<0
            commF(pop)=0;
        end
        if threadjN>numPy-1
            threadjN=0;
        elseif threadjN<0
            threadjN=numPy-1;
        end
        if threadkN>numPz-1
            commF(pop)=0;
        elseif threadkN<0
            commF(pop)=0;
        end
        
        if commF(pop)==1
            sendRank(pop)=threadiN+threadjN*(numPx)+threadkN*(numPy*numPx);
        end
    end
end

if periodicZ_==1 && periodicY_==1 && periodicX_==0
    for pop = 1:26
        threadiN=threadi+cx(pop);
        threadjN=threadj+cy(pop);
        threadkN=threadk+cz(pop);
        commF(pop)=1;
        if threadiN>numPx-1
            commF(pop)=0;
        elseif threadiN<0
            commF(pop)=0;
        end
        if threadjN>numPy-1
            threadjN=0;
        elseif threadjN<0
            threadjN=numPy-1;
        end
        if threadkN>numPz-1
            threadkN=0;
        elseif threadkN<0
            threadkN=numPz-1;
        end
        
        if commF(pop)==1
            sendRank(pop)=threadiN+threadjN*(numPx)+threadkN*(numPy*numPx);
        end
    end
end
if periodicZ_==1 && periodicY_==1 && periodicX_==1
    for pop = 1:26
        threadiN=threadi+cx(pop);
        threadjN=threadj+cy(pop);
        threadkN=threadk+cz(pop);
        commF(pop)=1;
        if threadiN>numPx-1
            threadiN=0;
        elseif threadiN<0
            threadiN=numPx-1;
        end
        if threadjN>numPy-1
            threadjN=0;
        elseif threadjN<0
            threadjN=numPy-1;
        end
        if threadkN>numPz-1
            threadkN=0;
        elseif threadkN<0
            threadkN=numPz-1;
        end
        
        if commF(pop)==1
            sendRank(pop)=threadiN+threadjN*(numPx)+threadkN*(numPy*numPx);
        end
    end
end

end

