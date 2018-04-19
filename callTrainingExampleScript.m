clear;clc;
for iijjkk=1:1
    clear gpsMeasStore imuMeasStore leverTrue Pstore ewxv0store ewxvPrevStore
    storemax=1e2;
    imuMeasMega=cell(100);
    gpsMeasMega=cell(100);
    leverMega=cell(100);
    stateMega=cell(100);
    for storeIndex=1:storemax
        tic;
        imu_training_examples
        ttm=toc;
        fprintf('%f%% complete. Run took %fs\n',storeIndex/storemax*100,ttm)
        imuMeasMega{storeIndex}=imuMeasStore;
        gpsMeasMega{storeIndex}=gpsMeasStore;
        leverMega{storeIndex}=Lab_true;
        stateMega{storeIndex}=statehist;
    end
    if iijjkk==1
        save truthData01.mat
        %save trainData13.mat
    elseif iijjkk==2
        save trainData02.mat
    elseif iijjkk==3
        save trainData03.mat
    elseif iijjkk==4
        save trainData04.mat
    elseif iijjkk==5
        save trainData05.mat
    elseif iijjkk==6
        save trainData06.mat
    elseif iijjkk==7
        save trainData07.mat
    elseif iijjkk==8
        save trainData08.mat
    elseif iijjkk==9
        save trainData09.mat
    elseif iijjkk==10
        save trainData10.mat
    elseif iijjkk==11
        save trainData11.mat
    elseif iijjkk==12
        save trainData12.mat
    elseif iijjkk==13
        save trainData13.mat
    elseif iijjkk==14
        save trainData14.mat
    elseif iijjkk==15
        save trainData15.mat
    elseif iijjkk==16
        save trainData16.mat
    elseif iijjkk==17
        save trainData17.mat
    elseif iijjkk==18
        save trainData18.mat
    elseif iijjkk==19
        save trainData19.mat
    elseif iijjkk==20
        save trainData20.mat
    elseif iijjkk==21
        save trainData21.mat
    elseif iijjkk==22
        save trainData22.mat
    elseif iijjkk==23
        save trainData23.mat
    elseif iijjkk==24
        save trainData24.mat
    elseif iijjkk==25
        save trainData25.mat
    elseif iijjkk==26
        save trainData26.mat
    elseif iijjkk==27
        save trainData27.mat
    elseif iijjkk==28
        save trainData28.mat
    elseif iijjkk==29
        save trainData29.mat
    elseif iijjkk==30
        save trainData30.mat
    end
end