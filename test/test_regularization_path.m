function test_suite = test_regularization_path
	buildFunctionHandleTestSuite(localfunctions);

end

function test_hcp25_participant
   
    datadir='/Volumes/MACBOOKUSB/Datasets/HCP_PTN1200/HCP_PTN1200_recon2';
    datafile='HCP1200_nodets_25ICA_200subjects.mat';
    hcpdata = load(fullfile(datadir,datafile));

    
	results = {};
    
    nsessions = 4;
    subject_idx = 1;
    for sessionno=1:nsessions
        sample_idx = (sessionno-1)*1200 + [1:1200];
        X = hcpdata.results.data(sample_idx,:,subject_idx);
        clear GGMobj
    	GGMobj = GGM(X,1,0);
    	[GGMobj results.sparsemle] = GGMobj.sparseMLE();
        figure(1)
        plot(results.sparsemle.Lrange,...
                results.sparsemle.sparsity,...
                '-',...
                'LineWidth',2,....
                'color',...
                [.8-.2*sessionno/nsessions .5 .2+.2*sessionno/nsessions]...
                );
        hold on;
        
        figure(2)
        subplot(nsessions/2,nsessions/2,sessionno);
        imagesc(GGMobj.ThetaPath(:,:,45));
    end
    
end

function test_acpi_participant
   
    datadir='/Volumes/MACBOOKUSB/Datasets/ACPI';
    datafile='dataset_acpi_nyu_task-rest_bold.mat';
    studydata = load(fullfile(datadir,datafile));

    
	results = {};
    
    nsessions = 2;
    sessionno = 1;
    subject_idx = 1;
    atlasno=2;
    for sessionno=1:2
        sample_idx = (sessionno-1)*180+[1:180];
        X = studydata.results{atlasno}.dataset{...
                                subject_idx}.Xdata(:,sample_idx)';
        X = standardize.successive_normalize(X);     
        
        ggm_options = estimator.create_options();
        ggm_options.nlambdas = 30;
        ggm_options.refit = false;
        p = size(X,2);
        results = estimator.fit(X,ggm_options);  
        sparsity = squeeze(sum(sum(...
                            abs(results.inverse_covariance_estimate)~=0)) ...
                            -p)/2;
        figure(1);
        plot(results.fit_options.path,...
                sparsity,...
                '-',...
                'LineWidth',2....
                );
        hold on;
        figure(2)
        subplot(1,nsessions,sessionno);
        imagesc(GGMobj.ThetaPath(:,:,20));
    end
    
end