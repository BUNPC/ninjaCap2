function testPlotPanels( o )


for iO = 1:length(o)
    
    v = o(iO).v;
    vOut = o(iO).vOut;
    e = o(iO).e;
    eOut = o(iO).eOut;
    
    if ~isfield(o(iO),'eSeamX')
        eSeamX = [];
    else
        eSeamX = o(iO).eSeamX;
        eSeamXpts = o(iO).eSeamXpts;
        eSeamXptsRot = o(iO).eSeamXptsRot;
    end

    figure(10+iO)
    plot(v(:,1),v(:,2),'k*', vOut(:,1),vOut(:,2),'ro')
    %set(gca,'ydir',"reverse")
    
    hold on
    for iE = 1:size(e,1)
        hl = plot(v(e(iE,:),1), v(e(iE,:),2), 'k-' );
        set(hl,'linewidth',2)
    end
    hold off
    
    hold on
    for iE = 1:size(eOut,1)
        plot(vOut(eOut(iE,:),1), vOut(eOut(iE,:),2), 'r:' );
    end
    hold off
    
    hold on
    plot( o(iO).xoutline, o(iO).youtline, 'r.-' )
    hold off
    
    hold on
    if iO==1 & isfield(o(2),'eSeamX')
        % left side
        for iS = 1:length(o(2).eSeamX)
            if o(2).eSeamX(iS)==1
                plot(o(2).eSeamXptsRot(iS,[1 3]), o(2).eSeamXptsRot(iS,[2 4]),'b.-')
                if o(2).eSeamXptsRot(iS,5)~=0 | o(2).eSeamXptsRot(iS,6)~=0
                    plot(o(2).eSeamXptsRot(iS,[1 5]), o(2).eSeamXptsRot(iS,[2 6]),'c.-')
                end
            end
        end
        
        % right side
        for iS = 1:length(o(3).eSeamX)
            if o(3).eSeamX(iS)==1
                plot(o(3).eSeamXptsRot(iS,[1 3]), o(3).eSeamXptsRot(iS,[2 4]),'b.-')
                if o(3).eSeamXptsRot(iS,5)~=0 | o(3).eSeamXptsRot(iS,6)~=0
                    plot(o(3).eSeamXptsRot(iS,[1 5]), o(3).eSeamXptsRot(iS,[2 6]),'c.-')
                end
            end
        end
        
    else
        
        for iS = 1:length(eSeamX)
            if eSeamX(iS)==1
                
                plot(eSeamXpts(iS,[1 3]), eSeamXpts(iS,[2 4]),'b.-')
                if eSeamXpts(iS,5)~=0 | eSeamXpts(iS,6)~=0
                    plot(eSeamXpts(iS,[1 5]), eSeamXpts(iS,[2 6]),'c.-')
                end
                
            end
        end
    end
    hold off
    
end

