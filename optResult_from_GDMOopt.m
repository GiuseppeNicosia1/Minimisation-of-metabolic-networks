function resObj = optResult_from_GDMOopt(GDMOobj,epsilon)

resObj = optResult(GDMOobj.results_folder, GDMOobj.gen, GDMOobj.optMode, epsilon, GDMOobj.V, GDMOobj.M, GDMOobj.plotFlag, GDMOobj.flagCorner, GDMOobj.utopianBio,GDMOobj.utopianSynth);

end
