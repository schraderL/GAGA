# Analyzing default repeat annotations

## Introduction
The basic TE annotation pipeline for GAGA looks as follows:

From Zijun:
>We run TE pipeline with the following steps:
> 1) Identification of known TEs
>We first identified known TEs using RepeatMasker (version 4.1.0) against the RepBase TE library (version RepBase21.12), and then executed RepeatProteinMask, which identifies TEs by aligning the genome sequence to the TE protein database in RepBase TE library.
> 2) De novo repeat prediction
>We constructed a de novo repeat library using RepeatModeler (version 2.0.1) with the default parameters. The generated results were consensus sequences and classification information for each repeat family. Then RepeatMasker (version 4.1.0) was run on the genome sequences, using the RepeatModeler consensus sequence as the library.  We also used ltr_finder (version v1.07) to identify long terminal repeat (LTR) retrotransposons in the genome.
> 3) Tandem repeats
>We predicted tandem repeats using TRF (version 4.09) with parameters set to "Match=2,Mismatch=7, Delta=7, PM=80, PI=10, Minscore=50, and MaxPeriod=2000".


Thus, repeats are annotated twice, once with the RepBase TE library and once with the de novo repeat library. Hence, there is not a single RepeatMasker run to combine both. On the upside, the RepBase TE annotations have no strong bias due to differences in de novo repeats.



## Data storage
RepeatMasker outputs and de novo libraries are stored on ERDA (`GAGA/Gene_annotation/GAGA_annotations/RepeatMasker_output/`). This folder was downloaded to IEBs `global/homes`:

```bash
cd /global/homes/jg/schradel/data/GAGA/Repeats
lftp io.erda.dk -p 21 -e "mirror GAGA/Gene_annotation/GAGA_annotations/RepeatMasker_output/; bye"
md5sum RepeatMasker_output/*.tar.gz
```


```
3de4df5e9359c1d73f7f86a4fb456e48  GAGA-0001_repeatMasker_out.tar.gz
722e0dd90b832d3c2eb6d2193d64127e  GAGA-0003_repeatMasker_out.tar.gz
609e91b3e82819db2231478bfe3a2720  GAGA-0004_repeatMasker_out.tar.gz
170e90703acad7bebfa8c502d70d1eda  GAGA-0014_repeatMasker_out.tar.gz
42b9c103e04676308889c499ff3cef53  GAGA-0020_repeatMasker_out.tar.gz
9f2b0a9e1d7f93bf5234ba9fe13ecc49  GAGA-0024_repeatMasker_out.tar.gz
e90dba0d83a72143f4eb33d8a0540908  GAGA-0025_repeatMasker_out.tar.gz
599ee4644ce7aa66379a8b2660e8452c  GAGA-0026_repeatMasker_out.tar.gz
d361db6fb225cc08dd8e4b9bbd521007  GAGA-0028_repeatMasker_out.tar.gz
45e94d8b8575f03f07fa25937437d539  GAGA-0063_repeatMasker_out.tar.gz
a43d897199e92db73307441a18f8cc13  GAGA-0074_repeatMasker_out.tar.gz
c2fc61dadaa1a4f0c2cdeac7ff10b5ac  GAGA-0080_repeatMasker_out.tar.gz
25f2208e85c7866b88cfd8bf4afd5b14  GAGA-0082_repeatMasker_out.tar.gz
9afe4d1a0753473a931740ec1dc07e09  GAGA-0083_repeatMasker_out.tar.gz
31e37515927a24bee2063c0c6a1ba91b  GAGA-0084_repeatMasker_out.tar.gz
48be3eeece9cda43431605b45396facc  GAGA-0085_repeatMasker_out.tar.gz
d11beab07ed3091f533c8bab9fbfd848  GAGA-0087_repeatMasker_out.tar.gz
55d4e8f3a9b7e39354b2fc38b242ad84  GAGA-0090_repeatMasker_out.tar.gz
074e9392ab43db8bb970733a69279d9e  GAGA-0098_repeatMasker_out.tar.gz
b2c5e29b4da2f9a03c08086be6c794e5  GAGA-0099_repeatMasker_out.tar.gz
a65ba32736a4304e28d64c6dbd8a6d2f  GAGA-0103_repeatMasker_out.tar.gz
09066aced2c57bc1b1869613c58e9952  GAGA-0109_repeatMasker_out.tar.gz
9cedfd27bcdc0b1644a435de8d4c7d77  GAGA-0114_repeatMasker_out.tar.gz
f6db0e495d64c242dfd5babff2ed8801  GAGA-0165_new_repeatMasker_out.tar.gz
dfab0ccdbd2145534e98b91dfae23194  GAGA-0177_repeatMasker_out.tar.gz
206396269e537e04d67248d5cb3e5724  GAGA-0187_repeatMasker_out.tar.gz
24abc9b1a54e0ed15607e991d237069a  GAGA-0198_repeatMasker_out.tar.gz
25e9d65e3a89c6eb2fd2e896655b7095  GAGA-0199_repeatMasker_out.tar.gz
9a9101a044754b3d401961b770a994fc  GAGA-0200_repeatMasker_out.tar.gz
ab40636b27dcc06077a4717f20a16929  GAGA-0221_repeatMasker_out.tar.gz
39cb37d8a6fb1fb292e667fcdaf1e6ff  GAGA-0222_repeatMasker_out.tar.gz
55bc99ec078f13ce6648deefc7dd367e  GAGA-0223_repeatMasker_out.tar.gz
f63b95364f1569d079a57b0d997a3747  GAGA-0224_repeatMasker_out.tar.gz
2daa6b4ae05528d729ed472a00ef069d  GAGA-0229_repeatMasker_out.tar.gz
4bf12bcfe637ee0676571f2c411d1668  GAGA-0245_repeatMasker_out.tar.gz
769761ab8693090057350e69c3ca1b3c  GAGA-0246_repeatMasker_out.tar.gz
27587b0d32f74cc7db8d73d7e3c7cce8  GAGA-0256_repeatMasker_out.tar.gz
b3f8bae4097f5ae61abaa17f78c4368a  GAGA-0266_repeatMasker_out.tar.gz
aa534fe5acefbcd9e043162917de751f  GAGA-0275_repeatMasker_out.tar.gz
eb61e8a63dd011fa28eacda57d39ac41  GAGA-0288_repeatMasker_out.tar.gz
12dd88f3f018d1725ad21c2866454bf0  GAGA-0300_repeatMasker_out.tar.gz
dc442ef6b996c6e5239cb6b8731f577c  GAGA-0301_repeatMasker_out.tar.gz
2a023e95a6cf5eb2f2a9b4ca7e26ad62  GAGA-0302_repeatMasker_out.tar.gz
c768ddea1d24b1a21b498e2e2533d1bd  GAGA-0303_repeatMasker_out.tar.gz
78bc72a916e4002ca7a5a06e589caac3  GAGA-0304_repeatMasker_out.tar.gz
f8c32efd1f450846577786b239ed52c9  GAGA-0306_repeatMasker_out.tar.gz
450519fd88e2ece4860057d3ec7da5df  GAGA-0307_repeatMasker_out.tar.gz
4f33f317d66431b216afbb89f1569bfa  GAGA-0328_repeatMasker_out.tar.gz
36639b55fe2cc28bf8cf9048bf74dde3  GAGA-0330_repeatMasker_out.tar.gz
ba8e8d619f6e838cfa8af5a8bf90db82  GAGA-0332_repeatMasker_out.tar.gz
84bebb2740a103efd4ff911dbc049089  GAGA-0333_repeatMasker_out.tar.gz
1d2d10a1504bcc16c04e637537c7f0e1  GAGA-0334_repeatMasker_out.tar.gz
883413a4197b579ba030bbee483464ea  GAGA-0335_repeatMasker_out.tar.gz
8542225553a472a020243a477cb9dab4  GAGA-0336_repeatMasker_out.tar.gz
293ef8b4b799109f8cafd8533ac63d31  GAGA-0337_repeatMasker_out.tar.gz
cf1a3174f882681fc50090333f8b4d24  GAGA-0338_repeatMasker_out.tar.gz
1a5b48eb559d015bb77121a81f1c5588  GAGA-0340_repeatMasker_out.tar.gz
bc15abddd3df86e5c3c88ee3007164fc  GAGA-0341_repeatMasker_out.tar.gz
96e4502e8854c21169aa7c45783211e4  GAGA-0343_repeatMasker_out.tar.gz
83ba289d7307a34a8e1ac39d4e65eb25  GAGA-0346_repeatMasker_out.tar.gz
661dcfb55789dea55fb4b6b8b5351800  GAGA-0350_repeatMasker_out.tar.gz
e6a3bb86d8e74c56e3dd11bada8c0930  GAGA-0351_repeatMasker_out.tar.gz
8e78efeb3eb5a1ffdbb0e9de6ae0d249  GAGA-0352_repeatMasker_out.tar.gz
0e876755942f327722d2f78cb9286beb  GAGA-0353_repeatMasker_out.tar.gz
ef614b05ed991181001c1ff1bfb49f7e  GAGA-0354_repeatMasker_out.tar.gz
36287603d5f1cc53d0b2462ace4f0bea  GAGA-0356_repeatMasker_out.tar.gz
926e22c078f184e0ca9aa5ab6b29069e  GAGA-0358_repeatMasker_out.tar.gz
c9a74831762bc43ee1811bdfb051007e  GAGA-0359_repeatMasker_out.tar.gz
f9d70c811aa44cbacda133e503f55bad  GAGA-0360_repeatMasker_out.tar.gz
4e96afec75a37afa13a19dc3932594de  GAGA-0361_repeatMasker_out.tar.gz
81a49374304769bd84338fc8d98dcc66  GAGA-0362_repeatMasker_out.tar.gz
c4164bd0a1856224743a7a1c9dc97861  GAGA-0363_repeatMasker_out.tar.gz
cf6e31364b4efe6cbadfe30eebe67e30  GAGA-0364_repeatMasker_out.tar.gz
e37cfdaa0aece98853134fd272e364d1  GAGA-0365_repeatMasker_out.tar.gz
d49ee34956a3ed0c5e6843f523a5b41b  GAGA-0374_repeatMasker_out.tar.gz
e6073a9d6d4c65b00aa2c36b4987d52f  GAGA-0376_repeatMasker_out.tar.gz
d14d72067fe4aab14fd5b24f02159fff  GAGA-0378_repeatMasker_out.tar.gz
73e0f9c373ea1e4a3221397fe530661e  GAGA-0379_repeatMasker_out.tar.gz
510013c604e7423bc814d7e735daa0a8  GAGA-0380_repeatMasker_out.tar.gz
dfa755bf4f2e3ace66e58e5e1de75a0d  GAGA-0382_repeatMasker_out.tar.gz
d7cbc2b592d626799cd28e9162cee67f  GAGA-0384_repeatMasker_out.tar.gz
accdce2b59ff35c10a37efc734549477  GAGA-0391_repeatMasker_out.tar.gz
8b43d6e6ef09d1373fdd6eafd8e52025  GAGA-0393_repeatMasker_out.tar.gz
952f04f13f7e5783144f24e645c1cb06  GAGA-0395_repeatMasker_out.tar.gz
35732e185cfb47b977da1e327cc89c15  GAGA-0396_repeatMasker_out.tar.gz
4f881c8a81aedff4c5783ab6fa226cfa  GAGA-0401_repeatMasker_out.tar.gz
e2bbd1c28415b58a65b512d6563cde11  GAGA-0404_repeatMasker_out.tar.gz
d29c9abb1a68840ed6e2743d419639f8  GAGA-0405_repeatMasker_out.tar.gz
759988b3511ab970aeda9ddaead5c41d  GAGA-0406_repeatMasker_out.tar.gz
f3e2cbd37b06a32d858b06f94197fedb  GAGA-0407_repeatMasker_out.tar.gz
31528c26dd11bf2945c223dab49f05bc  GAGA-0408_repeatMasker_out.tar.gz
56fe3590e5a75ea9d43ebfce56d0f7b8  GAGA-0413_repeatMasker_out.tar.gz
75bd479fe2dce141532515579026463b  GAGA-0454_repeatMasker_out.tar.gz
8c0c71836e1b4449f46710d846a5fc10  GAGA-0463_repeatMasker_out.tar.gz
d1830c0786110686b54c526208b946de  GAGA-0485_repeatMasker_out.tar.gz
6f29b9d259aad43ac41815db38d87810  GAGA-0491_repeatMasker_out.tar.gz
cfeb66e5c2a87c7b3ed0af62278d0ad6  GAGA-0494_repeatMasker_out.tar.gz
109588005e1720e58b2b9ac2ecdf0ae9  GAGA-0495_repeatMasker_out.tar.gz
f076dd10c7ab09681b03110cf23d01e6  GAGA-0502_repeatMasker_out.tar.gz
3510ed5f49a601ecc5510067cd973062  GAGA-0503_repeatMasker_out.tar.gz
472b149131a7c2e98cf10d47d51cd413  GAGA-0505_repeatMasker_out.tar.gz
9adc93283bd6fe19e8881f766fb58cc4  GAGA-0510_repeatMasker_out.tar.gz
6da93be3424f5653bdcdbd0bfb8b7177  GAGA-0511_repeatMasker_out.tar.gz
67f09c3f60a70202bb75d457552a1ea9  GAGA-0512_repeatMasker_out.tar.gz
831635154d25b142210548d5170123ea  GAGA-0513_repeatMasker_out.tar.gz
6f86be772ffcc8ca0086d100ba183e55  GAGA-0515_repeatMasker_out.tar.gz
d40b84d9a517c882f85898bd74b5abb4  GAGA-0517_repeatMasker_out.tar.gz
d57503481f0e200e7c2a6ba4e8cf4d7a  GAGA-0520_repeatMasker_out.tar.gz
0f83c741f40a464c382a01defe3840d4  GAGA-0521_repeatMasker_out.tar.gz
72e0093f484398a3fd5003694d443192  GAGA-0522_repeatMasker_out.tar.gz
66f36e62922750ce515ccbeb3bcb8d76  GAGA-0524_repeatMasker_out.tar.gz
9ab10025e257582f2e3aedbd19b76432  GAGA-0527_repeatMasker_out.tar.gz
6108e8233f09cdb601fb298d1d867b92  GAGA-0528_repeatMasker_out.tar.gz
89e0c6a391c64764998484ed36ebcf20  GAGA-0530_repeatMasker_out.tar.gz
3c9e1d04b892175431059df7750f6614  GAGA-0531_repeatMasker_out.tar.gz
9eb5191b56fb255171f78854cc6584a4  GAGA-0532_repeatMasker_out.tar.gz
c0df0c3342f28b2f1ae78528b796d35b  GAGA-0533_repeatMasker_out.tar.gz
14897af54725339139388f324fbaa25b  GAGA-0534_repeatMasker_out.tar.gz
ae0d054bab070aff0daade5a191a9520  GAGA-0535_repeatMasker_out.tar.gz
d32493c9bf06a66fecf7b1f22993a28e  GAGA-0536_repeatMasker_out.tar.gz
705b5cd6af47adc60f4f5194976054dd  GAGA-0537_repeatMasker_out.tar.gz
a72c3074c03766a63bf5f2175b883a31  GAGA-0538_repeatMasker_out.tar.gz
3df80632efb66bb68bf47cd11d84b6bf  GAGA-0539_repeatMasker_out.tar.gz
eba2a39ee0baf74cfb9bcb28c0956af8  GAGA-0540_repeatMasker_out.tar.gz
9a2eb724d5cb52c4c196f0cf0fac4912  GAGA-0541_repeatMasker_out.tar.gz
70ac77fef5f72c823dd5c12eb03af36c  GAGA-0543_repeatMasker_out.tar.gz
4d381b16d1ae911e25900fa1bc0f0ce0  GAGA-0550_repeatMasker_out.tar.gz
21e5b8089e71bc654e0642f8b15f99d6  GAGA-0552_repeatMasker_out.tar.gz
2b485c038bb4df5b4d48f69c995971e0  GAGA-0553_repeatMasker_out.tar.gz
bc1de9b4bf390cb7d932c57bc73e59b0  GAGA-0554_repeatMasker_out.tar.gz
0a4636bb8e16982a39aa04feaebde746  GAGA-0577_repeatMasker_out.tar.gz
b0d3fc657ee9419b26021b5a08dcae3c  GAGA-0578_repeatMasker_out.tar.gz
94667c1166b748ac3b3704e00b475bc5  GAGA-0579_repeatMasker_out.tar.gz
21b647bae3af54732d5aaaee11e7b445  GAGA-0580_repeatMasker_out.tar.gz
5f48e9979de350b2e5cd1720654a02dd  NCBI-0001_repeatMasker_out.tar.gz
2e26b30daf7ca2af46a63c3eb465dbfb  NCBI-0002_repeatMasker_out.tar.gz
358a97b72b78fae290d1fd8af5841e60  NCBI-0003_repeatMasker_out.tar.gz
15416125c8cc187cdc723bd9acf42ea5  NCBI-0004_repeatMasker_out.tar.gz
747a897a8a8fc1c0c7815160b792d9a5  NCBI-0005_repeatMasker_out.tar.gz
741ff9a368c0f559b3fd65be61e49f7e  NCBI-0006_repeatMasker_out.tar.gz
a8172105626fb8b110b7806b4073eb48  NCBI-0007_repeatMasker_out.tar.gz
9dceb3daeb3392d36063924eb2671f19  NCBI-0008_repeatMasker_out.tar.gz
4337ebbf430d808a64aff0864f31a599  NCBI-0009_repeatMasker_out.tar.gz
820ddeb86574627549fce9d42bc42376  NCBI-0010_repeatMasker_out.tar.gz
1bc12fd58caadb0952973363211c1fa5  NCBI-0011_repeatMasker_out.tar.gz
d4710d10f603bd6b7dff99586e490e47  NCBI-0012_repeatMasker_out.tar.gz
864b8c80133e99d36ec8bc459c741ea7  NCBI-0013_repeatMasker_out.tar.gz
790246ff2b070fb9e1c09507fa60f00e  NCBI-0014_repeatMasker_out.tar.gz
75d18ab06946208edb9f4955b39563d8  NCBI-0015_repeatMasker_out.tar.gz
9bb40ad7726b06b6985e9ae87698c9b5  NCBI-0016_repeatMasker_out.tar.gz
5a87ca31fc43bfd77b6615870e090c6e  NCBI-0017_repeatMasker_out.tar.gz
e7a30f9348d870e10fb04448eadc659a  OUT-0001_repeatMasker_out.tar.gz
de2565038a999beec084d566e54a7f40  OUT-0002_repeatMasker_out.tar.gz
```
