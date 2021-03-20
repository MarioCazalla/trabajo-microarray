#Obtenemos el numero de lineas del fichero targets.txt
line_number=$(wc -l targets.txt | awk '{print $1}')

#Le restamos 1 para obtener el numero de muestras
sample_number=$(($line_number-1))

#Creamos el fichero .cls y añadimos el número de muestras en la primera linea
echo $sample_number > ./GSEA_files/GSE18198.cls

#Contamos cuantas clases hay
classes_number=$(awk 'NR!=1{print $3}' targets.txt | sort | uniq | wc -l)
#Añadimos un espacio, el número de clases, otro espacio y un 1 tras el número de muestras
sed "s/$/ $classes_number 1/" ./GSEA_files/GSE18198.cls > ./GSEA_files/GSE18198_V2.cls
#Borramos el archivo anterior y renombramos el nuevo
rm ./GSEA_files/GSE18198.cls; mv ./GSEA_files/GSE18198_V2.cls ./GSEA_files/GSE18198.cls

#Añadimos la almohadilla de la segunda línea
echo "#" >> ./GSEA_files/GSE18198.cls

#Guardamos los nombres de las clases separados por espacios
classes_names=$(awk '$1!="FileName"{print $3}' targets.txt | sort | uniq)
#Susitutimos los saltos de linea por espacios
classes_names=$(echo $classes_names | sed "s/\n/ /")
echo $classes_names
#Sustituimos el asterisco por la linea completa con asterisco + clases
awk -v var="$classes_names" 'NR==2 {gsub("#","# "var)}; {print}' ./GSEA_files/GSE18198.cls > ./GSEA_files/GSE18198_V2.cls
rm ./GSEA_files/GSE18198.cls; mv ./GSEA_files/GSE18198_V2.cls ./GSEA_files/GSE18198.cls

touch cls_third_line.txt
for line in targets.txt
do
temp=$(echo $line | awk '{print $3}')
if [ "$temp" == "DMSO" ]
then 
    echo 0 >> cls_third_line.txt
elif [ "$temp" == "SAHM1" ]
then
    echo 1 >> cls_third_line.txt
fi
done
