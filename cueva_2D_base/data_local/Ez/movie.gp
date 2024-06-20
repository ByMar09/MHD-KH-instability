# Definir el nombre base y el rango de números
filename = "2d_data_ez"
min_ext = 1
max_ext = 30

# Configurar el formato de salida de la gráfica pm3d
set terminal pngcairo

# Establecemos el fondo negro
set object 1 rectangle from screen 0, 0 to screen 1, 1 behind fillcolor rgb "black" fillstyle solid noborder

# Establecemos el titulo
set title "Resistive Rotor" font "Times-Roman,22"  tc rgb "white"

# Establecemos la división de la malla
set xlabel "x" font "Times-Roman,18" tc rgb "white"
set xrange [ -0.500 : 0.500 ] noreverse nowriteback
set xtics  -0.5, 0.2, 0.5 font "Times-Roman,10"
set ylabel "y" rotate by 0  font "Times-Roman,18" tc rgb "white" 
set yrange [0.0 : 1] noreverse nowriteback
set ytics   0, 0.2, 1 font "Times-Roman,10"
set zlabel "Ez" 
set zrange [ -1: 1 ] noreverse nowriteback
set style line 50 lt 1 lc rgb "white" lw 2
set border ls 50
set cblabel "" font "Times-Roman,15" #tc rgb "white"
set cbrange [-1:1]
set cbtics  -1,0.5,1

# Configurar el estilo pm3d
set pm3d map

# Configurar el rango de colores
set palette rgb 33,13,10

# Configurar la grilla
set grid

# Definir la función para graficar un archivo
plot_file(ext) = sprintf('< tar xfO %s.tar %s_%04d.dat.bz2 | bunzip2 -c', filename, filename, ext) 

# Graficar todos los archivos en el rango y guardar cada imagen
do for [ext=min_ext:max_ext] {
    set output sprintf("grafica_%04d.png", ext)
    splot plot_file(ext)
}

# Convertir las imágenes en un archivo .avi utilizando ffmpeg
system("ffmpeg -y -r 10 -i grafica_%04d.png -vcodec mpeg4 -q:v 2 movie_ez.avi")
