# To run (replace "1-03" with the desired release):
# ./update_maracluster.sh 1-03
release_long=$1
release=${release_long:0:4}
echo $release_long
echo $release
# Update the linux binary
wget https://github.com/statisticalbiotechnology/maracluster/releases/download/rel-${release_long}/tarball.tar.gz

tar -xzvf tarball.tar.gz
tar -xzvf maracluster-v${release}-linux-amd64.tar.gz

cp maracluster-v${release}-linux-amd64/bin/maracluster simsi_transfer/utils/maracluster/linux64/

rm tarball.tar.gz
rm maracluster-v${release}-linux-amd64.tar.gz

rm -rf maracluster-v${release}-linux-amd64


# Update the windows binary
wget https://github.com/statisticalbiotechnology/maracluster/releases/download/rel-${release_long}/win64.zip
unzip win64.zip

file-roller -h maracluster-v${release}-win32-msvc-x86.exe

cp maracluster-v${release}-win32-msvc-x86/bin/maracluster.exe simsi_transfer/utils/maracluster/win64/

rm win64.zip
rm maracluster-v${release}-win32-msvc-x86.exe
rm -rf maracluster-v${release}-win32-msvc-x86
