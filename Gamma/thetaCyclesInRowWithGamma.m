function [slowInRowHist, fastInRowHist] = thetaCyclesInRowWithGamma(indexSlow, indexFast, thetaIndices)




N = length(thetaIndices);

theta_cycles_slow_fast_indices = zeros(N,2);

theta_cycles_slow_fast_indices(indexSlow,1) = 1;
theta_cycles_slow_fast_indices(indexFast,2) = 1;

count1 = 0;
count2 = 0;
count3 = 0;
count4 = 0;
count5 = 0;
count6 = 0;
count7 = 0;
count8 = 0;
count9 = 0;
count10 = 0;

for k = 1:N-10
    if theta_cycles_slow_fast_indices(k,1) == 1 
        if theta_cycles_slow_fast_indices(k+1,1) ~= 1
            count1 = count1 + 1;
        end
        
        if theta_cycles_slow_fast_indices(k+1,1) == 1
            if theta_cycles_slow_fast_indices(k+2,1) ~= 1
                count2 = count2 + 1;
            end
            
            if theta_cycles_slow_fast_indices(k+2,1) == 1
                if theta_cycles_slow_fast_indices(k+3,1) ~= 1
                    count3 = count3 + 1;
                end
                
                if theta_cycles_slow_fast_indices(k+3,1) == 1
                    if theta_cycles_slow_fast_indices(k+4,1) ~= 1
                        count4 = count4 + 1;
                    end
                    
                    if theta_cycles_slow_fast_indices(k+4,1) == 1
                        if theta_cycles_slow_fast_indices(k+5,1) ~= 1
                            count5 = count5 + 1;
                        end
                        
                        if theta_cycles_slow_fast_indices(k+5,1) == 1
                            if theta_cycles_slow_fast_indices(k+6,1) ~= 1
                                count6 = count6 + 1;
                            end
                            
                            if theta_cycles_slow_fast_indices(k+6,1) == 1
                                if theta_cycles_slow_fast_indices(k+7,1) ~= 1
                                    count7 = count7 + 1;
                                end
                                
                                if theta_cycles_slow_fast_indices(k+7,1) == 1
                                    if theta_cycles_slow_fast_indices(k+8,1) ~= 1
                                        count8 = count8 + 1;
                                    end
                                    
                                    if theta_cycles_slow_fast_indices(k+8,1) == 1
                                        if theta_cycles_slow_fast_indices(k+9,1) ~= 1
                                            count9 = count9 + 1;
                                        end
                                        
                                        if theta_cycles_slow_fast_indices(k+9,1) == 1
                                            if theta_cycles_slow_fast_indices(k+10,1) ~= 1
                                                count10 = count10 + 1;
                                            end
                                        end
                                    end
                                end
                            end
                        end    
                    end
                end
            end
        end
    end
end

slowInRowHist =[count1,count2,count3,count4,count5,count6,count7,count8,count9,count10];

slowInRowHist = slowInRowHist / (N-10);




count1 = 0;
count2 = 0;
count3 = 0;
count4 = 0;
count5 = 0;
count6 = 0;
count7 = 0;
count8 = 0;
count9 = 0;
count10 = 0;

for k = 1:N-10
    if theta_cycles_slow_fast_indices(k,2) == 1 
        if theta_cycles_slow_fast_indices(k+1,2) ~= 1
            count1 = count1 + 1;
        end
        
        if theta_cycles_slow_fast_indices(k+1,2) == 1
            if theta_cycles_slow_fast_indices(k+2,2) ~= 1
                count2 = count2 + 1;
            end
            
            if theta_cycles_slow_fast_indices(k+2,2) == 1
                if theta_cycles_slow_fast_indices(k+3,2) ~= 1
                    count3 = count3 + 1;
                end
                
                if theta_cycles_slow_fast_indices(k+3,2) == 1
                    if theta_cycles_slow_fast_indices(k+4,2) ~= 1
                        count4 = count4 + 1;
                    end
                    
                    if theta_cycles_slow_fast_indices(k+4,2) == 1
                        if theta_cycles_slow_fast_indices(k+5,2) ~= 1
                            count5 = count5 + 1;
                        end
                        
                        if theta_cycles_slow_fast_indices(k+5,2) == 1
                            if theta_cycles_slow_fast_indices(k+6,2) ~= 1
                                count6 = count6 + 1;
                            end
                            
                            if theta_cycles_slow_fast_indices(k+6,2) == 1
                                if theta_cycles_slow_fast_indices(k+7,2) ~= 1
                                    count7 = count7 + 1;
                                end
                                
                                if theta_cycles_slow_fast_indices(k+7,2) == 1
                                    if theta_cycles_slow_fast_indices(k+8,2) ~= 1
                                        count8 = count8 + 1;
                                    end
                                    
                                    if theta_cycles_slow_fast_indices(k+8,2) == 1
                                        if theta_cycles_slow_fast_indices(k+9,2) ~= 1
                                            count9 = count9 + 1;
                                        end
                                        
                                        if theta_cycles_slow_fast_indices(k+9,2) == 1
                                            if theta_cycles_slow_fast_indices(k+10,2) ~= 1
                                                count10 = count10 + 1;
                                            end
                                        end
                                    end
                                end
                            end
                        end    
                    end
                end
            end
        end
    end
end

fastInRowHist =[count1,count2,count3,count4,count5,count6,count7,count8,count9,count10];

fastInRowHist = fastInRowHist / (N-10);



