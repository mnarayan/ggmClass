function test_suite = test_stars
	buildFunctionHandleTestSuite(localfunctions);

end

function test_identity_graph
   
    p = 15;
    nlambdas = 10; 
    nresamples = 50; 
    grphs = {};
    for resampleno=1:nresamples
        grphs{resampleno,1} = repmat(eye(p),[p p nlambdas]); 
    end
    [scores best_lambda] = estimator.stars(grphs);
    assertVectorsAlmostEqual(scores,zeros(1,nlambdas),'Scores are not empty');
    assertEqual(best_lambda,1,'Fails Identity Test');
    
end


function test_dense_graph
   
    p = 15;
    nlambdas = 10; 
    nresamples = 50; 
    grphs = {};
    for resampleno=1:nresamples
        grphs{resampleno,1} = cat(3,repmat(eye(p),[p p nlambdas/2]),...
                                    repmat(ones(p,p),[p p nlambdas/2]) ...
                                    );
    end
    [scores best_lambda] = estimator.stars(grphs);    
    assertVectorsAlmostEqual(scores,zeros(1,nlambdas),'Scores are not empty');
    assertEqual(best_lambda,1,'Fails Dense Test');
end