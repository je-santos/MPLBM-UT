import pyvista.examples as examples
import pyvista as pv
import numpy as np
import os
import imageio
import re
import time
import sys


def natural_sort(l):
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
    return sorted(l, key=alphanum_key)


print("Creating GIF...")
cwd = os.getcwd()
images = []
filenames = os.listdir(cwd + '/gif_images/')
filenames = natural_sort(filenames)

for filename in filenames:
    images.append(imageio.imread('gif_images/' + filename))
imageio.mimsave('lbm_animation.gif', images)

print("Done!")
