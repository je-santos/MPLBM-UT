import os
import imageio
import re
import moviepy.editor as mp
import mplbm_utils as mplbm


def create_gif(anim_dir, save_name):

    print("Creating GIF...This may take a few minutes. Please be patient!")
    images = []
    filenames = os.listdir(anim_dir)
    filenames = mplbm.natural_sort(filenames)

    for filename in filenames:
        images.append(imageio.imread(anim_dir + filename))
    imageio.mimsave(anim_dir + "../" + save_name + "_lbm_animation.gif", images)

    print("Done!")


def create_mp4(anim_dir, save_name, speed_factor):

    print("Creating MP4...")
    gif = mp.VideoFileClip(anim_dir + "../" + save_name + "_lbm_animation.gif")
    gif = gif.speedx(factor=speed_factor)
    gif.write_videofile(anim_dir + "../" + save_name + "_lbm_animation.mp4")

    print("Done!")
