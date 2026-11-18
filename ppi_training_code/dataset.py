import torch
from torch.utils.data import Dataset, DataLoader
import numpy as np
import pickle
import random
from os.path import join, isfile

class Dataset_PPI_Interact(Dataset):
    # predict for all samples of all casp15 targets in each epoch validation, instead of only randomly select one sample for each target
    def __init__(self,
                 target_list,
                 esm_dir="/mnt/rna01/zhy/tools/PPLM-PPI_data/D_script/PPLM_feats/",
                 esm_flag=False,
                 esm_chain_dir="/mnt/rna01/junl/PPLM/PPLM_PPI/D_script/PPLM_ESM2_t33_chain/",
                 debug=False,
                 n_debug_target=100,
                 crop_size=512,
                 n_layer=12,
                 validation=False,
                 testing=False,
                 ):

        self.esm_dir = esm_dir
        self.esm_flag = esm_flag
        self.esm_chain_dir = esm_chain_dir
        self.testing = testing
        # self.feat_mode = feat_mode


        print("self.esm_dir:", self.esm_dir)

        self.n_layer = n_layer
        self.crop_size = crop_size

        self.samples = []
        for line in open(target_list, 'r').readlines():
            sample = line.strip().split()[0]

            if testing:
                # print(self.esm_dir, sample)
                esm_path = join(self.esm_dir, sample + ".pkl")
                if isfile(esm_path):
                    self.samples.append(sample)
            else:
                self.samples.append(sample)

            if debug and len(self.samples) >= n_debug_target:
                break

        if not validation or not testing:
            np.random.shuffle(self.samples)

        print("Dataset_PPI_Interact | number of targets:", len(self.samples))

    def __len__(self):
        return len(self.samples)

    def __getitem__(self, idx, pindex=-1):
        if torch.is_tensor(idx):
            idx = idx.tolist()

        sample = self.samples[idx]
        unsuccess = True
        while unsuccess:
            try:
                esm_path = join(self.esm_dir, sample + ".pkl")
                with open(esm_path, 'rb') as fr:
                    esm_data = pickle.load(fr)

                mean_embed_A = esm_data["mean_embed_A"]
                mean_embed_B = esm_data["mean_embed_B"]
                mean_attn_AA = esm_data["mean_attn_AA"]
                mean_attn_BB = esm_data["mean_attn_BB"]
                mean_inter_attn = esm_data["mean_inter_attn"]

                max_embed_A = esm_data["max_embed_A"]
                max_embed_B = esm_data["max_embed_B"]
                max_attn_AA = esm_data["max_attn_AA"]
                max_attn_BB = esm_data["max_attn_BB"]
                max_inter_attn = esm_data["max_inter_attn"]
                #
                min_embed_A = esm_data["min_embed_A"]
                min_embed_B = esm_data["min_embed_B"]
                min_attn_AA = esm_data["min_attn_AA"]
                min_attn_BB = esm_data["min_attn_BB"]
                min_inter_attn = esm_data["min_inter_attn"]

                # lenA = esm_data["lens"][0]
                # lenB = esm_data["lens"][1]

                label = esm_data["label"]

                if self.esm_flag:
                    chainA_name = sample.split("@")[0]
                    chainB_name = sample.split("@")[1]
                    esm_chainA_path = join(self.esm_chain_dir, chainA_name + ".fa.pkl")
                    esm_chainB_path = join(self.esm_chain_dir, chainB_name + ".fa.pkl")
                    with open(esm_chainA_path, 'rb') as fr:
                        esm_chainA_data = pickle.load(fr)
                    with open(esm_chainB_path, 'rb') as fr:
                        esm_chainB_data = pickle.load(fr)

                    chain_mean_embed_A = esm_chainA_data["mean_embed"]
                    chain_mean_embed_B = esm_chainB_data["mean_embed"]
                    chain_mean_attn_A = esm_chainA_data["mean_attn"]
                    chain_mean_attn_B = esm_chainB_data["mean_attn"]

                    chain_max_embed_A = esm_chainA_data["max_embed"]
                    chain_max_embed_B = esm_chainB_data["max_embed"]
                    chain_max_attn_A = esm_chainA_data["max_attn"]
                    chain_max_attn_B = esm_chainB_data["max_attn"]

                    chain_min_embed_A = esm_chainA_data["min_embed"]
                    chain_min_embed_B = esm_chainB_data["min_embed"]
                    chain_min_attn_A = esm_chainA_data["min_attn"]
                    chain_min_attn_B = esm_chainB_data["min_attn"]


                p = random.random()
                if p < 0.5:
                    data = {
                        # "sample": sample,
                        "mean_embed_A": mean_embed_A,
                        "mean_embed_B": mean_embed_B,
                        "mean_attn_AA": mean_attn_AA,
                        "mean_attn_BB": mean_attn_BB,
                        "mean_inter_attn": mean_inter_attn,
                        "max_embed_A": max_embed_A,
                        "max_embed_B": max_embed_B,
                        "max_attn_AA": max_attn_AA,
                        "max_attn_BB": max_attn_BB,
                        "max_inter_attn": max_inter_attn,
                        "min_embed_A": min_embed_A,
                        "min_embed_B": min_embed_B,
                        "min_attn_AA": min_attn_AA,
                        "min_attn_BB": min_attn_BB,
                        "min_inter_attn": min_inter_attn,
                        "label": label
                    }
                    if self.esm_flag:
                        data['chain_mean_embed_A'] = chain_mean_embed_A
                        data['chain_mean_embed_B'] = chain_mean_embed_B
                        data['chain_mean_attn_A'] = chain_mean_attn_A
                        data['chain_mean_attn_B'] = chain_mean_attn_B
                        data['chain_max_embed_A'] = chain_max_embed_A
                        data['chain_max_embed_B'] = chain_max_embed_B
                        data['chain_max_attn_A'] = chain_max_attn_A
                        data['chain_max_attn_B'] = chain_max_attn_B
                        data['chain_min_embed_A'] = chain_min_embed_A
                        data['chain_min_embed_B'] = chain_min_embed_B
                        data['chain_min_attn_A'] = chain_min_attn_A
                        data['chain_min_attn_B'] = chain_min_attn_B
                else:
                    data = {
                        # "sample": sample,
                        "mean_embed_A": mean_embed_B,
                        "mean_embed_B": mean_embed_A,
                        "mean_attn_AA": mean_attn_BB,
                        "mean_attn_BB": mean_attn_AA,
                        "mean_inter_attn": mean_inter_attn,
                        "max_embed_A": max_embed_B,
                        "max_embed_B": max_embed_A,
                        "max_attn_AA": max_attn_BB,
                        "max_attn_BB": max_attn_AA,
                        "max_inter_attn": max_inter_attn,
                        "min_embed_A": min_embed_B,
                        "min_embed_B": min_embed_A,
                        "min_attn_AA": min_attn_BB,
                        "min_attn_BB": min_attn_AA,
                        "min_inter_attn": min_inter_attn,
                        "label": label
                    }
                    if self.esm_flag:
                        data['chain_mean_embed_A'] = chain_mean_embed_B
                        data['chain_mean_embed_B'] = chain_mean_embed_A
                        data['chain_mean_attn_A'] = chain_mean_attn_B
                        data['chain_mean_attn_B'] = chain_mean_attn_A
                        data['chain_max_embed_A'] = chain_max_embed_B
                        data['chain_max_embed_B'] = chain_max_embed_A
                        data['chain_max_attn_A'] = chain_max_attn_B
                        data['chain_max_attn_B'] = chain_max_attn_A
                        data['chain_min_embed_A'] = chain_min_embed_B
                        data['chain_min_embed_B'] = chain_min_embed_A
                        data['chain_min_attn_A'] = chain_min_attn_B
                        data['chain_min_attn_B'] = chain_min_attn_A
                unsuccess = False


            except Exception as e:
                print("Warning: ", sample, esm_path)
                new_idx = np.random.choice(np.arange(len(self.samples)))
                sample = self.samples[new_idx]

        return data
