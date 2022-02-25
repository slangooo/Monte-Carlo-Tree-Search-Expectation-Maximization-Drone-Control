x = [mvnrnd([5 10],reshape([12 8],1,2).*eye(2),100); mvnrnd([6 9],reshape([4 8],1,2).*eye(2),100);mvnrnd([2 2],reshape([4 8],1,2).*eye(2),100)];
m=15;
sigma_in=1.*500.*ones(1,m);
softmaxScale= 4;
maxUsersPerCluster=40;
maxIter=100;
sinr_threshold=2;
maxNumOfServedUsers = 0;
for i1=500:-100:0
    for i2=500:-100:0
        for i3=500:-100:0
            for i4=500:-100:0
                for i5=500:-100:0
                    for i6=500:-100:0
                        for i7=500:-100:0
                            for i8=500:-100:0
                                for i9=500:-100:0
                                    for i10=500:-100:0
                                        for i11=500:-100:0
                                            for i12=500:-100:0
                                                for i13=500:-100:0
                                                    for i14=500:-100:0
                                                        for (i15=0:5)
                                                            sigma_in=[i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12,i13,i14,i15*100];
                                                            [ efficiencyScore,allocationScore,reliabilityScore, maxSINR, numberOfServedUsers(i15+1),iter] = EM_PULL_CNST_SIGMA1...
                                                                (x,m,sigma_in,softmaxScale,maxUsersPerCluster,maxIter,sinr_threshold);
                                                        end
                                                        [maxNumberPerIteration,idx] = max(numberOfServedUsers);
                                                        sigma_in(end) = (idx-1)*100;
                                                        if maxNumberPerIteration > maxNumOfServedUsers
                                                            maxSigma = sigma_in;
                                                            maxNumOfServedUsers = maxNumberPerIteration;
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
        end
    end
end



