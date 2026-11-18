import sys
import argparse
import os

network_path = os.path.dirname(__file__)
# print("   PPLM | PPI-Contact: ", network_path)
sys.path.insert(0, network_path)
from PPI_model import *
from dataset import Dataset_PPI_Interact

os.environ['MKL_THREADING_LAYER'] = 'GNU'
os.environ["CUDA_VISIBLE_DEVICES"] = "0,1,2,3,4,5,6,7"

import numpy as np
import torch
from torch.utils.data import Dataset, DataLoader

def main():
    parser = argparse.ArgumentParser(description="Error predictor network trainer",
                                     epilog="v0.0.1")

    parser.add_argument("out_path_prefix",
                        action="store",
                        help="Location of folder to save checkpoints to.")

    parser.add_argument("--mean_path",
                        "-mean_path",
                        type=str,
                        default='',
                        help="top5 model pkl list path",
                        )

    parser.add_argument("--max_path",
                        "-max_path",
                        type=str,
                        default='',
                        help="top5 model pkl list path",
                        )

    parser.add_argument("--min_path",
                        "-min_path",
                        type=str,
                        default='',
                        help="top5 model pkl list path",
                        )

    parser.add_argument("--order_mode",
                        "-order_mode",
                        type=str,
                        default="both",
                        help="order_mode: random, normal, reverse, both",
                        )

    parser.add_argument("--e_45",
                        "-e_45",
                        action="store_true",
                        default=False,
                        help="pplm_e_45 model (Default: False)")

    parser.add_argument("--linear_all",
                        "-linear_all",
                        action="store_true",
                        default=False,
                        help="linear_all (Default: False)")

    parser.add_argument("--mean_max_min_update",
                        "-mean_max_min_update",
                        action="store_true",
                        default=False,
                        help="mean_max_min_update (Default: False)")

    parser.add_argument("--only_inter_attn",
                        "-only_inter_attn",
                        action="store_true",
                        default=False,
                        help="only_inter_attn (Default: False)")

    parser.add_argument("--embed",
                        "-embed",
                        action="store_true",
                        default=False,
                        help="embed (Default: False)")

    parser.add_argument("--gpu_id",
                        "-gpu",
                        type=int,
                        default=0,
                        help="gpu device specified",
                        )

    parser.add_argument("--debug",
                        "-debug",
                        action="store_true",
                        default=False,
                        help="Debug mode (Default: False)")


    args = parser.parse_args()

    ################################################## Define model ##################################################
    print("Defining network model...")

    assigned_device = "cuda:" + str(args.gpu_id)
    device = assigned_device if torch.cuda.is_available() else "cpu"
    print("   device: ", device)

    esm_flag = False

    model_mean_paths = []
    if args.mean_path != '':
        for line in open(args.mean_path, 'r').readlines():
            model_mean_paths.append(line.strip())

    model_max_paths = []
    if args.max_path != '':
        for line in open(args.max_path, 'r').readlines():
            model_max_paths.append(line.strip())

    model_min_paths = []
    if args.min_path != '':
        for line in open(args.min_path, 'r').readlines():
            model_min_paths.append(line.strip())

    print("Mean models:", model_mean_paths)
    print("Max models:", model_max_paths)
    print("Min models:", model_min_paths)

    model_mean_list = []
    for model_path in model_mean_paths:
        model = PPI_inter_intra_attn_embed_single_pooling()
        checkpoint = torch.load(model_path, map_location=device)
        model.load_state_dict(checkpoint["model_state_dict"])
        model.to(device)
        model.eval()
        total_num = sum(p.numel() for p in model.parameters())
        trainable_num = sum(p.numel() for p in model.parameters() if p.requires_grad)
        print("Number of Network Parameters | Total:", total_num, "Trainable", trainable_num)

        model_mean_list.append(model)

    model_max_list = []
    for model_path in model_max_paths:
        model = PPI_inter_intra_attn_embed_single_pooling()
        checkpoint = torch.load(model_path, map_location=device)
        model.load_state_dict(checkpoint["model_state_dict"])
        model.to(device)
        model.eval()
        total_num = sum(p.numel() for p in model.parameters())
        trainable_num = sum(p.numel() for p in model.parameters() if p.requires_grad)
        print("Number of Network Parameters | Total:", total_num, "Trainable", trainable_num)

        model_max_list.append(model)

    model_min_list = []
    for model_path in model_min_paths:
        model = PPI_inter_intra_attn_embed_single_pooling()
        checkpoint = torch.load(model_path, map_location=device)
        model.load_state_dict(checkpoint["model_state_dict"])
        model.to(device)
        model.eval()
        total_num = sum(p.numel() for p in model.parameters())
        trainable_num = sum(p.numel() for p in model.parameters() if p.requires_grad)
        print("Number of Network Parameters | Total:", total_num, "Trainable", trainable_num)

        model_min_list.append(model)


    test_sets = ['mouse_test', 'fly_test', 'worm_test', 'yeast_test', 'ecoli_test']

    for test_set in test_sets:

        testset_path = "/mnt/rna01/zhy/tools/PPLM-PPI_data/D_script/test_set." + test_set + ".txt"
        print("testset_path:", testset_path)

        test_samples = Dataset_PPI_Interact(testset_path, esm_flag=esm_flag, order_mode=args.order_mode, debug=args.debug, n_debug_target=1000, e_45=args.e_45)
        test_dataloader = DataLoader(test_samples, batch_size=512, num_workers=8, prefetch_factor=4, pin_memory=True)  # , prefetch_factor=8)
        ##################################################################################################################

        #############################################################################################
        #################################### Validation stage #########################################
        #############################################################################################
        # print("Validation...")
        # model.eval()  # turn off training mode
        temp_loss = {"total": []}
        with torch.no_grad():  # wihout tracking gradients
            total_preds_ave = torch.Tensor()
            total_labels = torch.Tensor()
            total_preds_ave = total_preds_ave.to(device)
            total_labels = total_labels.to(device)

            # output = [["name", "label", "prediction"]]
            out_statistic_path = args.out_path_prefix + ".preds.statistic"
            out_stat = open(out_statistic_path, 'w')
            out_stat.write(f"output_statistic: {out_statistic_path}")

            output_path = args.out_path_prefix + "." + test_set + ".preds.csv"
            output = open(output_path, 'w')
            output.write("name,label,prediction\n")
            print("output_path:", output_path)
            out_stat.write(f"output_csv: {output_path}")


            for j, data in enumerate(test_dataloader):
                labels = data['label'].to(device)
                sample = data['sample']

                predictions_list = []

                # if args.mode == 'mean_max' or args.mode == 'mean_min' or args.mode == 'mean_max_min':
                for model in model_mean_list:
                    mean_embed_A = data['mean_embed_A'].to(device)
                    mean_embed_B = data['mean_embed_B'].to(device)
                    mean_attn_AA = data['mean_attn_AA'].to(device)
                    mean_attn_BB = data['mean_attn_BB'].to(device)
                    mean_inter_attn = data['mean_inter_attn'].to(device)

                    predictions = model(mean_inter_attn, mean_attn_AA, mean_attn_BB, mean_embed_A, mean_embed_B)
                    if args.order_mode == 'both':
                        predictions2 = model(mean_inter_attn, mean_attn_BB, mean_attn_AA, mean_embed_B, mean_embed_A)
                        predictions = (predictions + predictions2) / 2

                    predictions_list.append(predictions)

                # if args.mode == 'mean_max' or args.mode == 'max_min' or args.mode == 'mean_max_min':
                for model in model_max_list:
                    max_embed_A = data['max_embed_A'].to(device)
                    max_embed_B = data['max_embed_B'].to(device)
                    max_attn_AA = data['max_attn_AA'].to(device)
                    max_attn_BB = data['max_attn_BB'].to(device)
                    max_inter_attn = data['max_inter_attn'].to(device)

                    predictions = model(max_inter_attn, max_attn_AA, max_attn_BB, max_embed_A, max_embed_B)
                    if args.order_mode == 'both':
                        predictions2 = model(max_inter_attn, max_attn_BB, max_attn_AA, max_embed_B, max_embed_A)
                        predictions = (predictions + predictions2) / 2

                    predictions_list.append(predictions)

                # if args.mode == 'mean_min' or args.mode == 'max_min' or args.mode == 'mean_max_min':
                for model in model_min_list:
                    min_embed_A = data['min_embed_A'].to(device)
                    min_embed_B = data['min_embed_B'].to(device)
                    min_attn_AA = data['min_attn_AA'].to(device)
                    min_attn_BB = data['min_attn_BB'].to(device)
                    min_inter_attn = data['min_inter_attn'].to(device)

                    predictions = model(min_inter_attn, min_attn_AA, min_attn_BB, min_embed_A, min_embed_B)
                    if args.order_mode == 'both':
                        predictions2 = model(min_inter_attn, min_attn_BB, min_attn_AA, min_embed_B, min_embed_A)
                        predictions = (predictions + predictions2) / 2

                    predictions_list.append(predictions)


                predictions = torch.stack(predictions_list)
                predictions = torch.mean(predictions, dim=0)

                # print("predictions:", predictions.shape, labels.shape, len(sample))
                # output.append([sample, labels.detach().cpu().item(), predictions.detach().cpu().item()])
                for idx in range(len(sample)):
                    name = sample[idx]
                    label = labels[idx].detach().cpu().item()
                    pred = predictions[idx].detach().cpu().item()
                    output.write(name + "," + str(label) + "," + str(pred) + "\n")

                loss = torch.nn.functional.binary_cross_entropy(predictions.squeeze(), labels.float())

                temp_loss["total"].append(loss.item())
                ###########################################################################################

                total_preds_ave = torch.cat((total_preds_ave, predictions.squeeze()), 0)
                total_labels = torch.cat((total_labels, labels), 0)

            total_preds_ave = total_preds_ave.cpu().numpy().flatten()
            total_labels = total_labels.cpu().numpy().flatten()

            output.close()
            # with open(output_path, "w", newline="", encoding="utf-8") as f:
            #     writer = csv.writer(f)
            #     writer.writerow(output)

            precision, recall, accuracy, F1_score, Top10, Top20, Top50, AUC_ROC, AUC_PR = evaluate(total_labels, total_preds_ave)

            out_stat.write(f"==========: {test_set} | Loss: {np.mean(temp_loss['total'])} | Accuracy: {precision} {recall} {accuracy} {F1_score} {AUC_ROC} {AUC_PR}")
            print("==========", test_set, "| Loss:", np.mean(temp_loss["total"]), "| Accuarcy:", precision, recall, accuracy, F1_score, Top10, Top20, Top50, AUC_ROC, AUC_PR)


if __name__ == "__main__":
    main()