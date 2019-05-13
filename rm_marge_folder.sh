bart_dir='/sfs/bart'
for user_key in $bart_dir/usercase/*;
do
    if [ -d $user_key ];
    then
        find $user_key/marge* -maxdepth 0 -type d -ctime +10 -exec rm -rf {} + >> $bart_dir/log/clean_dir.log 2>&1; # remove those marge_* files
    fi
done

