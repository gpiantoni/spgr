
for subject in {'EM09','MG17','MG33','MG37','MG61','MG63'}
do link2dir /home/gio/recordings/${subject}/mri/proc/freesurfer/ /home/gio/projects/spgr/group/fsaverage/${subject}
done

export SUBJECTS_DIR='/home/gio/projects/spgr/group/fsaverage'

ln -s /opt/freesurfer/subjects/fsaverage_sym ${SUBJECTS_DIR}/fsaverage_sym

for subject in {'EM09','MG17','MG33','MG37','MG61','MG63'}
do xhemireg --s $subject
done

for subject in {'EM09','MG17','MG33','MG37','MG61','MG63'};
do surfreg --s $subject --t fsaverage_sym --lhrh;
surfreg --s $subject --t fsaverage_sym --xhemi --lhrh;
done

make_average_subject --out myatlas.i1 \
  --surf-reg fsaverage_sym.sphere.reg \
  --subjects EM09 MG17 MG33 MG37 MG61 MG63 \
  --xhemi \
  --no-vol --template-only 

for subject in {'EM09','MG17','MG33','MG37','MG61','MG63'}
do  surfreg --s $subject --t myatlas.i1 --lhrh; 
surfreg --s $subject --t myatlas.i1 --xhemi --lhrh;
done

make_average_subject --out myatlas.i2 \
  --surf-reg myatlas.i1.sphere.reg \
  --subjects EM09 MG17 MG33 MG37 MG61 MG63 \
  --xhemi \
  --no-vol --template-only 

for subject in {'EM09','MG17','MG33','MG37','MG61','MG63'}
do surfreg --s $subject --t myatlas.i2 --lhrh;
surfreg --s $subject --t myatlas.i2 --xhemi --lhrh;
done

make_average_subject --out myatlas.i3 \
  --surf-reg myatlas.i2.sphere.reg \
  --subjects EM09 MG17 MG33 MG37 MG61 MG63 \
  --xhemi


