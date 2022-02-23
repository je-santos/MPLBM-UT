import os
import imageio
import re
import moviepy.editor as mp
import sys
sys.path.append('../../python_utils/')
from parse_input_file import *


def natural_sort(l):
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
    return sorted(l, key=alphanum_key)

def create_gif(anim_dir):

    print("Creating GIF...This may take a few minutes. Please be patient!")
    images = []
    filenames = os.listdir(anim_dir)
    filenames = natural_sort(filenames)

    for filename in filenames:
        images.append(imageio.imread(anim_dir + filename))
    imageio.mimsave(anim_dir + "../" + save_name + "_lbm_animation.gif", images)

    print("Done!")


def create_mp4(anim_dir, speed_factor):

    print("Creating MP4...")
    gif = mp.VideoFileClip(anim_dir + "../" + save_name + "_lbm_animation.gif")
    gif = gif.speedx(factor=speed_factor)
    gif.write_videofile(anim_dir + "../" + save_name + "_lbm_animation.mp4")

    print("Done!")


input_file = 'input.yml'
inputs = parse_input_file(input_file)  # Parse inputs
inputs['input output']['simulation directory'] = os.getcwd()  # Store current working directory

anim_dir = inputs['input output']['output folder'] + 'animation/'
save_name = inputs['domain']['geom name']
create_gif(anim_dir)
create_mp4(anim_dir, speed_factor=1)  # speed_factor = 1 means no slow down or speed up
