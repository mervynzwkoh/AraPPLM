import sys
import argparse
import os

network_path = os.path.dirname(__file__)
# print("   PPLM | PPI-Contact: ", network_path)
sys.path.insert(0, network_path)

os.environ['MKL_THREADING_LAYER'] = 'GNU'
os.environ["CUDA_VISIBLE_DEVICES"] = "0,1,2,3,4,5,6,7"

from os.path import isfile, isdir, join
import numpy as np
import torch
from time import *
import random
from torch.utils.data import Dataset, DataLoader

from PPI_model import *
from dataset import Dataset_PPI_Interact

def main():
    parser = argparse.ArgumentParser(description="Error predictor network trainer",
                                     epilog="v0.0.1")

    parser.add_argument("folder",
                        action="store",
                        help="Location of folder to save checkpoints to.")

    parser.add_argument("--cross_round",
                        "-cross_round",
                        type=int,
                        default=0,
                        help="cross validation round",
                        )

    parser.add_argument("--seed_round",
                        "-seed_round",
                        type=int,
                        default=0,
                        help="gpu device specified",
                        )

    parser.add_argument("--linear_all",
                        "-linear_all",
                        action="store_true",
                        default=False,
                        help="linear_all (Default: False)")

    parser.add_argument("--local_rank",
                        type=int,
                        help="# of epochs (path over all proteins) to train for (Default: 200)")

    parser.add_argument("--mode",
                        "-mode",
                        type=str,
                        default="PPI_inter_intra_attn_embed_mean_max",
                        help="feature mode: "
                             "PPI_inter_intra_attn_embed_mean_max_min,"
                             "PPI_inter_intra_attn_embed_mean_max,"
                             "PPI_inter_intra_attn_embed_mean,"
                             "PPI_inter_intra_attn_embed_max,"
                             "PPI_inter_intra_attn_embed_min,"
                             "PPI_inter_intra_attn_mean_max,"
                             "PPI_intra_attn_mean_max,"
                             "PPI_inter_attn_mean_max,"
                             "PPI_embed_mean_max,"
                             "PPI_intra_attn_embed_mean_max"
                             "PPI_inter_intra_attn_mean,"
                             "PPI_intra_attn_mean,"
                             "PPI_inter_attn_mean,"
                             "PPI_embed_mean,"
                             "PPI_intra_attn_embed_mean,"
                             "PPLM_PPI_inter_intra_attn_embed_mean_max_min,"
                             "PPLM_PPI_inter_intra_attn_embed_mean_max,"
                        )

    parser.add_argument("--batch_size",
                        "-batch_size", action="store",
                        type=int,
                        default=512,
                        help="batch_size: 512, 128, 32")

    parser.add_argument("--learning_rate",
                        "-lr",
                        type=float,
                        default=0.00005,
                        help="gpu device specified",
                        )

    parser.add_argument("--seed_add",
                        "-seed",
                        type=int,
                        default=0,
                        help="gpu device specified",
                        )


    parser.add_argument("--gpu_id",
                        "-gpu",
                        type=int,
                        default=0,
                        help="gpu device specified",
                        )

    parser.add_argument("--n_layer",
                        "-n_layer",
                        type=int,
                        default=12,
                        help="gpu device specified",
                        )

    parser.add_argument("--epochs",
                        "-e", action="store",
                        type=int,
                        default=15,
                        help="# of epochs (path over all proteins) to train for (Default: 200)")

    parser.add_argument("--debug",
                        "-debug",
                        action="store_true",
                        default=False,
                        help="Debug mode (Default: False)")

    parser.add_argument("--silent",
                        "-s",
                        action="store_true",
                        default=False,
                        help="Run in silent mode (Default: False)")


    args = parser.parse_args()
    # ##################################################################################################################

    print("PPLM: PPLM_PPI | mode:", args.mode)

    seeds = [26240761, 73591024, 48136297, 95317406, 62780953, 58302741, 94716038, 12649570, 80421359, 39287564]
    # seed = 26240761 + args.seed_add
    seed = seeds[args.seed_round] + args.seed_add
    seed_everything(seed)

    print("random.seed(): ", random.seed())
    print("np.random.get_state()[1][0]: ", np.random.get_state()[1][0])
    print("torch.initial_seed(): ", torch.initial_seed())
    print("torch.cuda.initial_seed(): ", torch.cuda.initial_seed())
    # ##################################################################################################################

    cross_validation_round = str(args.cross_round)
    seed_round = str(args.seed_round)
    loss_name = args.folder + ".round_" + seed_round + ".cv_" + cross_validation_round #+ ".lr_" + str(args.learning_rate)
    if args.linear_all:
        loss_name = loss_name + ".linear_all"
    out_recording_file = loss_name + ".recording"

    trianset_path = "/mnt/rna01/zhy/tools/PPLM-PPI_data/D_script/human.tsv_correct.cross_validation_" + cross_validation_round + ".trainset"
    validset_path = "/mnt/rna01/zhy/tools/PPLM-PPI_data/D_script/human.tsv_correct.cross_validation_" + cross_validation_round + ".testset"

    ################################################## Define model ##################################################
    print("Defining network model...")

    assigned_device = "cuda:" + str(args.gpu_id)
    device = assigned_device if torch.cuda.is_available() else "cpu"
    print("   device: ", device)
    esm_flag = False

    # args.mode == "PPI_inter_intra_attn_embed_mean" or args.mode == "PPI_inter_intra_attn_embed_max" or args.mode == "PPI_inter_intra_attn_embed_min":
    model = PPI_inter_intra_attn_embed_single_pooling()

    model.to(device)

    total_num = sum(p.numel() for p in model.parameters())
    trainable_num = sum(p.numel() for p in model.parameters() if p.requires_grad)
    print("Number of Network Parameters | Total:", total_num, "Trainable", trainable_num)

    ### Creating folder ##################################################
    if isfile(join(loss_name, "model.pkl")):
        print("Loading a checkpoint...")
        checkpoint = torch.load(join(loss_name, "model.pkl"), map_location=device)
        model.load_state_dict(checkpoint["model_state_dict"])
        epoch = checkpoint["epoch"] + 1
        train_loss = checkpoint["train_loss"]
        valid_loss = checkpoint["valid_loss"]
        best_models = checkpoint["best_models"]

        if not args.silent: print("Restarting at epoch", epoch)
    else:
        print("Training a new model...")
        epoch = 1
        train_loss = {"total": []}
        valid_loss = {"total": []}
        best_models = []
        if not isdir(loss_name):
            if not args.silent: print("Creating a new dir at", loss_name)
            os.makedirs(loss_name)


    if epoch > args.epochs:
        exit()

    out_recording = open(out_recording_file, "a+")
    out_recording.write(loss_name + "\n")
    out_recording.write("PPLM: PPLM_PPI | mode:" + args.mode + "\n")
    out_recording.write("random.seed: " + str(26240761 + args.seed_add) + "\n")
    out_recording.write("trianset_path: " + trianset_path + "\n")
    out_recording.write("validset_path: " + validset_path + "\n")
    out_recording.write("Number of Network Parameters | Total:" + str(total_num) + " Trainable: " + str(trainable_num) + "\n")
    out_recording.write("learning rate: " + str(args.learning_rate) + "\n")
    out_recording.write("batch_size: " + str(args.batch_size) + "\n")



    ################################################## Loading data ##################################################
    if not args.silent: print("Loading dataset...")
    print("trian_pre_load_path:", trianset_path)
    print("valid_pre_load_path:", validset_path)

    train_samples = f(trianset_path, esm_flag=esm_flag, debug=args.debug, n_debug_target=5000)
    train_dataloader = DataLoader(train_samples, batch_size=args.batch_size, shuffle=True, num_workers=8, prefetch_factor=4, pin_memory=True)  # , prefetch_factor=8)
    valid_samples = Dataset_PPI_Interact(validset_path, esm_flag=esm_flag, validation=True, debug=args.debug, n_debug_target=1000)
    valid_dataloader = DataLoader(valid_samples, batch_size=args.batch_size, num_workers=8, prefetch_factor=4, pin_memory=True)  # , prefetch_factor=8)
    ##################################################################################################################

    best_AUC_PR = [-1, -1, -1, -1, -1]
    best_AUC_ROC = [-1, -1, -1, -1, -1]
    best_F1_score = [-1, -1, -1, -1, -1]

    ################################################## Network optimizer ##################################################print("Fixed learning rate...:", LR)
    optimizer = torch.optim.AdamW(model.parameters(), lr=args.learning_rate, weight_decay=0.01)  # 0.001
    start_epoch = epoch

    # 混合双精度训练
    # scaler = GradScaler()
    for epoch in range(start_epoch, args.epochs + 1):
        begin_time_0 = time()
        model.train(True)  # 开启网络模型训练模式
        temp_loss = {"total": []}
        # print('Training...')
        for i, data in enumerate(train_dataloader):
            #################################### Feature load #######################################
            optimizer.zero_grad()

            labels = data['label'].to(device)

            if args.mode == "PPI_inter_intra_attn_embed_mean":
                mean_embed_A = data['mean_embed_A'].to(device)
                mean_embed_B = data['mean_embed_B'].to(device)
                mean_attn_AA = data['mean_attn_AA'].to(device)
                mean_attn_BB = data['mean_attn_BB'].to(device)
                mean_inter_attn = data['mean_inter_attn'].to(device)

                predictions = model(mean_inter_attn, mean_attn_AA, mean_attn_BB, mean_embed_A, mean_embed_B)

            elif args.mode == "PPI_inter_intra_attn_embed_max":
                max_embed_A = data['max_embed_A'].to(device)
                max_embed_B = data['max_embed_B'].to(device)
                max_attn_AA = data['max_attn_AA'].to(device)
                max_attn_BB = data['max_attn_BB'].to(device)
                max_inter_attn = data['max_inter_attn'].to(device)

                predictions = model(max_inter_attn, max_attn_AA, max_attn_BB, max_embed_A, max_embed_B)

            elif args.mode == "PPI_inter_intra_attn_embed_min":
                min_embed_A = data['min_embed_A'].to(device)
                min_embed_B = data['min_embed_B'].to(device)
                min_attn_AA = data['min_attn_AA'].to(device)
                min_attn_BB = data['min_attn_BB'].to(device)
                min_inter_attn = data['min_inter_attn'].to(device)

                predictions = model(min_inter_attn, min_attn_AA, min_attn_BB, min_embed_A, min_embed_B)


            loss = torch.nn.functional.binary_cross_entropy(predictions.squeeze(), labels.float())


            loss.backward()
            optimizer.step()

            temp_loss["total"].append(loss.item())


            ###########################################################################################

            # if (i + 1) % 1000 == 0 or (i + 1) == len(train_dataloader):
            #     sys.stdout.write(
            #             "\rEpoch: %2d, [%2d/%2d] | Train_loss: %.3f | Time: %4d\n"
            #             % (epoch, i+1, len(train_dataloader), loss.item(), time() - begin_time_0)
            #         )

        # print("End training")

        train_loss["total"].append(np.mean(temp_loss["total"]))


        #############################################################################################
        #################################### Validation stage #########################################
        #############################################################################################
        # print("Validation...")
        temp_loss = {"total": []}
        model.eval()  # turn off training mode
        with torch.no_grad():  # wihout tracking gradients
            total_preds_ave = torch.Tensor()
            total_labels = torch.Tensor()
            total_preds_ave = total_preds_ave.to(device)
            total_labels = total_labels.to(device)

            for j, data in enumerate(valid_dataloader):
                labels = data['label'].to(device)

                if args.mode == "PPI_inter_intra_attn_embed_mean":
                    mean_embed_A = data['mean_embed_A'].to(device)
                    mean_embed_B = data['mean_embed_B'].to(device)
                    mean_attn_AA = data['mean_attn_AA'].to(device)
                    mean_attn_BB = data['mean_attn_BB'].to(device)
                    mean_inter_attn = data['mean_inter_attn'].to(device)

                    predictions = model(mean_inter_attn, mean_attn_AA, mean_attn_BB, mean_embed_A, mean_embed_B)

                elif args.mode == "PPI_inter_intra_attn_embed_max":
                    max_embed_A = data['max_embed_A'].to(device)
                    max_embed_B = data['max_embed_B'].to(device)
                    max_attn_AA = data['max_attn_AA'].to(device)
                    max_attn_BB = data['max_attn_BB'].to(device)
                    max_inter_attn = data['max_inter_attn'].to(device)

                    predictions = model(max_inter_attn, max_attn_AA, max_attn_BB, max_embed_A, max_embed_B)

                elif args.mode == "PPI_inter_intra_attn_embed_min":
                    min_embed_A = data['min_embed_A'].to(device)
                    min_embed_B = data['min_embed_B'].to(device)
                    min_attn_AA = data['min_attn_AA'].to(device)
                    min_attn_BB = data['min_attn_BB'].to(device)
                    min_inter_attn = data['min_inter_attn'].to(device)

                    predictions = model(min_inter_attn, min_attn_AA, min_attn_BB, min_embed_A, min_embed_B)

                loss = torch.nn.functional.binary_cross_entropy(predictions.squeeze(), labels.float())

                temp_loss["total"].append(loss.item())
                ###########################################################################################

                total_preds_ave = torch.cat((total_preds_ave, predictions.squeeze()), 0)
                total_labels = torch.cat((total_labels, labels), 0)

                # if (j + 1) % 1000 == 0 or (j + 1)== len(valid_dataloader):
                #     sys.stdout.write(
                #             "\rEpoch: %2d, [%2d/%2d] | Train_loss: %.3f | Time: %4d\n"
                #             % (epoch, j+1, len(valid_dataloader), loss.item(), time() - begin_time_0)
                #         )

            valid_loss["total"].append(np.mean(temp_loss["total"]))


            total_preds_ave = total_preds_ave.cpu().numpy().flatten()
            total_labels = total_labels.cpu().numpy().flatten()
            # print("XXX:", total_preds_ave.shape, total_labels.shape)

            # print("total_preds_ave:", total_preds_ave)
            # print("total_labels:", total_labels)

            precision, recall, accuracy, F1_score, Top10, Top20, Top50, AUC_ROC, AUC_PR = evaluate(total_labels, total_preds_ave)

            #################################### Saving models #######################################

        #############################################################################################
        ############################## Saving current epoch model ###################################
        #############################################################################################
        torch.save({
            'epoch': epoch,
            'model_state_dict': model.state_dict(),
            # 'optimizer_state_dict': optimizer.state_dict(),
            'train_loss': train_loss,
            'valid_loss': valid_loss,
            'best_models': best_models,
        }, join(loss_name, "model_" + str(epoch) + ".pkl"))

        torch.save({
            'epoch': epoch,
            'model_state_dict': model.state_dict(),
            # 'optimizer_state_dict': optimizer.state_dict(),
            'train_loss': train_loss,
            'valid_loss': valid_loss,
            'best_models': best_models,
        }, join(loss_name, "model.pkl"))


        loss_record = "========== Loss Epoch: " + str(epoch) + "/" + str(args.epochs) + " lr: " + str(optimizer.param_groups[0]['lr']) + " time: " + str(time() - begin_time_0) + " | Train_loss: " + str(train_loss["total"][-1]) + " Valid_loss: " + str(valid_loss["total"][-1])
        accu_record = "========== Accuarcy Epoch: " + str(epoch) + "/" + str(args.epochs) + " | " + str(precision) + " " + str(recall) + " " + str(accuracy) + " " + str(F1_score) + " " + str(Top10) + " " + str(Top20) + " " + str(Top50) + " " + str(AUC_ROC) + " " + str(AUC_PR)

        out_recording.write(loss_record + "\n")
        out_recording.write(accu_record + "\n")

        print(loss_record)
        print(accu_record)

    out_recording.close()


def seed_everything(seed=26240761):
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False
    torch.backends.cudnn.enabled = False

def worker_init_fn(worker_id, rank, seed=26240761):
    worker_seed = rank + seed
    random.seed(worker_seed)
    np.random.seed(worker_seed)
    torch.manual_seed(worker_seed)

if __name__ == "__main__":
    main()