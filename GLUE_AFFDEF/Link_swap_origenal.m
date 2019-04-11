x=[0.35 0.0230 0.417 0.01200 1.30000 20.18 0.0393 0.3110 0.00510 1.09000 8.62];
for jj=1:size(x,1)
    fid1 = fopen('Swap.swp' , 'r');
    fid2 = fopen('SwapCopy.swp' ,'w');

    for i=1:312
        s = fgets(fid1);
        fprintf(fid2 , '%s', s);
    end

    f1New = x(jj,1);

    f1Old = fscanf(fid1,' COFRED = %f ! Soil evaporation coefficient of Black, [0..1 cm/d1/2, R],');
    fprintf(fid2 , ' COFRED = %4.2f ! Soil evaporation coefficient of Black, [0..1 cm/d1/2, R],\n',  f1New);

    for i=1:39
        s = fgets(fid1);
        fprintf(fid2 , '%s', s);
    end

    f3Old= fscanf(fid1, '%f ', 8 );
    f3New = [1 2 3 4 5 6 7 8]; f3New = f3Old ;
    f3New(2)= x(jj,2); f3New(3)= x(jj,3); f3New(4)= x(jj,4); f3New(5)= x(jj,5); f3New(6)= x(jj,6);
    fprintf(fid2, '         %g  %7.4f%7.3f   %7.5f  %7.5f%7.2f    %7.5f  %7.5f\r\n', f3New(1), f3New(2),...
    f3New(3), f3New(4), f3New(5), f3New(6), f3New(7), f3New(8) );

    f4Old= fscanf(fid1, '%f', 8 );
    f4New = [2 4 6 8 10 12 14 16]; f4New = f4Old ;
    f4New(2)= x(jj,7); f4New(3)= x(jj,8); f4New(4)= x(jj,9); f4New(5)= x(jj,10); f4New(6)= x(jj,11);
    fprintf(fid2, '         %g  %7.4f%7.3f   %7.5f  %7.5f%7.2f    %7.5f  %7.5f', f4New(1), f4New(2),...
    f4New(3), f4New(4), f4New(5), f4New(6), f4New(7), f4New(8) );

    while(1)
        s = fgets(fid1);
        if(s<0)
            break;
        end
        fprintf(fid2 , '%s', s);
    end

    fclose(fid1);
    fclose(fid2);

    delete ('Swap.swp');
    a = 'SwapCopy.swp'; b = 'Swap.swp'; 
    eval(['!rename ' a ' ' b]);

    % Call model to generate simulated data
    dos SWAP;

    % read moisture output 
    fid = fopen('result.dat');
    r = textscan(fid, '%*s %*f %f %f %f %f %*f','HeaderLines', 8);
    fclose(fid);
    rr = cell2mat(r);
    rr(11,:)=[];
    simdata = zeros(40,1);
    simdata(:)=rr;
end