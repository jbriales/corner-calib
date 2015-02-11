classdef CExperimentState < handle
    properties
        per  = 0 % Current percentage
        counter = 0
        total
        t1 % First time
    end
    
    methods
        function this = CExperimentState( N )
            this.total = N;
            this.t1 = tic;
        end
        
        function update( this )
            this.counter = this.counter + 1;
            if floor(this.counter/this.total*100) > this.per
                this.per = floor(this.counter/this.total*100);
                fprintf('%d/100\t%d of %d\tEstimated remaining %d\n',...
                    this.per,this.counter,this.total,...
                    floor((this.total-this.counter)*toc(this.t1)/this.counter));
            end
        end
        
        function reset( this )
            this.per = 0;
            this.counter = 0;
        end
    end
end